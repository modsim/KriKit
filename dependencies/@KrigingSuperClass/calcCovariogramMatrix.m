function [ ] = calcCovariogramMatrix(obj)
% []=calcCovariogramMatrix()
%
% Calculates the extended (co)variance matrix and, if UseInverse=true, its
% inverse.
% The actual (co)variance matrix contains the evaluation of the
% (co)variogram model for each pair-wise combination of provided sample
% data. This Kriging-(Co)variogram-Matrix is symmetrically extended by
% evaluations of given basis functions(ignoring linear parameters).
% Concequently, the final matrix has the size of
% (nExperiments+nBasisFunctions)X(nExperiments+nBasisFunctions) 
%
% Before calcCovariogramMatrix() is called, nonlinear parameters should be
% calculated using "solveLeastSquareBasisFct"/"solveLeastSquareBasisFctGA".
%
% Note: If the (co)variogram model contains the parameter sigmaError, this
% parameter is only apllied on the diangonal of the matrix. This is done in
% order to "simulate" measurement noise
%
% You can set: -
% - CovariogramModelChoice ... define which covariogram is used
% - UseInverse ... Decide if the inverse of the (co)variance matrix should
%                  be calculated. The inverse is needed in "prediction()".
%                  If "UseInverse=false", Gaussian elemenation is
%                  alternatively used.
% - BasisFct ... set basis function which are used for later on for
%                prediction
%
% You can get:
% - CovariogramMatrix ... saves the covariance model values for
%                         pairwise combinations of the experimental data +
%                         evaluations of the basis functions
% - Variogram ... saves the variance model values for
%                 pairwise combinations of the experimental data +
%                 evaluations of the basis functions
% - InvVariogram ... inverse of "Variogram"
% - InvCovariogramMatrix ... inverse of "CovariogramMatrix"
% 
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

%% Initialization
% Check if Basis Funtions are defined
if (obj.nBasisFct<1)
    error('No basis function was defined. Use the defineBasisFunction(nBasisFct,type,varargin) for defining basis functions!')
end

% General Variogram Matrix
    % Initialization
obj.Variogram = zeros(obj.nExperiments+obj.nBasisFct,obj.nExperiments+obj.nBasisFct);


% All combinations
[a,b] = ndgrid(1:obj.nExperiments,1:obj.nExperiments);
combinationMatrix = [b(:),a(:)];

if obj.CovariogramModelChoice>0&&isempty(obj.DistInput)
    error('DistInput is empty please run obj.estimateVariance')
end

%% Calulate the distance and (Co)variogram matrix
if obj.CovariogramUsesEuclideanDistance
    delta = sqrt(sum((bsxfun(@minus,obj.InputData(combinationMatrix(:,1),:),obj.InputData(combinationMatrix(:,2),:)).^2),2));
    delta = reshape(delta,sqrt(size(delta,1)),sqrt(size(delta,1)));

    % Calculate the variogram values for theses distances
    obj.Variogram(1:obj.nExperiments,1:obj.nExperiments) =  obj.CovarModel(zeros(1,size(obj.DistInput,2)),1) - obj.CovarModel(delta,0);
elseif obj.CovariogramUsesAbsoluteDistance
    % Calculate the variogram values for the input distances (for
    % each input variable). 
    
    % "delta" size: 
    % (nInputDataXnInputData)X(nInputVar,1)X(nInputDataXnInputData)
    % and is interpretated as
    % [pairwise distance of with respect to input variable 1]X[scalar for input variable index to which is comapred with]X[pairwise distance of with respect to input variable 2]
    delta=abs(bsxfun(@minus,obj.InputData(combinationMatrix(:,1),:),obj.InputData(combinationMatrix(:,2),:)));
    nRowsColumns=sqrt(size(delta,1));
    d = ones(nRowsColumns,nRowsColumns);
    for iN = 1:nRowsColumns
        for iV = 1:obj.nInputVar
           d(:,(iN-1)*obj.nInputVar+iV)=delta(nRowsColumns*(iN-1)+1:nRowsColumns*(iN),iV);
        end
    end
    delta = reshape(d,nRowsColumns,obj.nInputVar,nRowsColumns);
    obj.Variogram(1:obj.nExperiments,1:obj.nExperiments) =  reshape(obj.CovarModel(zeros([1,size(obj.DistInput,2),1]),1) - obj.CovarModel(delta,0),nRowsColumns,nRowsColumns);
else
    error('Either"CovariogramUsesEuclideanDistance" or "CovariogramUsesAbsoluteDistance" have to be "true"')
end

    % Diagonal has to be 0
obj.Variogram = obj.Variogram - bsxfun(@times,diag(obj.Variogram),eye(size(obj.Variogram,1)));


% Addiitonal extend by ones
% Extend Matrix by Basis Functions
for iBasis = 1 : obj.nBasisFct
        basis = obj.BasisFct{iBasis}(obj.BasisFctParameters,obj.getInputData);
        obj.Variogram(1:obj.getnExperiments,obj.nExperiments+iBasis) = basis;
        obj.Variogram(obj.nExperiments+iBasis,1:obj.getnExperiments) = basis;
end

% Make sure that the last elements in the right corner are zeros
obj.Variogram(obj.nExperiments+1:obj.nExperiments+obj.nBasisFct,obj.nExperiments+1:obj.nExperiments+obj.nBasisFct) = zeros(obj.nBasisFct);

% Convert to CavoriogramMatrix
obj.CovariogramMatrix = obj.Variogram;
covariance_zero  = ones(obj.nExperiments,obj.nExperiments)*obj.CovarModel(zeros(1,size(obj.DistInput,2)),1);
obj.CovariogramMatrix(1:obj.nExperiments,1:obj.nExperiments) = covariance_zero-obj.Variogram(1:obj.nExperiments,1:obj.nExperiments);

%% Calculate Inverse if needed
% If the inverse is used for inv(A)*b instead of A\b
if obj.UseInverse
    % Calculate Inverse
    if obj.nBasisFctParameters == 0&&obj.getUseSimpleKriging==1
        obj.InvVariogram=inv(obj.Variogram(1:obj.nExperiments,1:obj.nExperiments));
        obj.InvCovariogramMatrix=inv(obj.CovariogramMatrix(1:obj.nExperiments,1:obj.nExperiments));
    else
        % Use Explicit matrix (If number of basis function is big, this can be beneficial)
        % Calculate the inverse directly
        obj.InvVariogram  = inv(obj.Variogram);
        obj.InvCovariogramMatrix = inv(obj.CovariogramMatrix);
    end
end

end

% =============================================================================
%  KriKit - Kriging toolKit
%  
%  Copyright 2014-2016: Lars Freier(1), Eric von Lieres(1)
%                                      
%    (1)Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
%  
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
