function [] = estimateVariance(obj)
% [] = estimateVariance()
% Calculates an estimation of the variance between measurements as
% dependency of there distance to each other. The calculation of the distance
% depends on the chosen model. So it maybe the euclidean norm (scalar) or 
% the absolute difference (vector) of the input data. The distance is saved 
% in the vector/Matrix DistInput and the estimation of the Variance in the
% DistSquareOutput
%
% You can set: -
% 
% You can get:
%- DistInput ... contains the distances of input variables to each other
%            if euclidean: dimension is (nExperiments*(nExperiments-1)/2)X1
%            if absolute: dimension is (nExperiments*(nExperiments-1)/2)XnInputVar
% - DistSquareOutput ... Matern estimation of the variance of the output.
%                    This is needed for the variogram fitting apporach
%                    (used for covariogram parameter estimation)
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.
%%
if obj.nExperiments<=1

    error('At least two measurements must be provided');
end

    if obj.CovariogramUsesAbsoluteDistance&&obj.CovariogramUsesEuclideanDistance
        error('"CovariogramUsesEuclideanDistance" and "CovariogramUsesAbsoluteDistance" are "true"')
    end
        
    if obj.CovariogramUsesEuclideanDistance
        obj.DistInput = [];
        obj.VarianceEstimation = [];

        % -------------------------------------------------------------
        % Calculate for every combination of the experiments. Since
        % here this distances are only used for variogramm the data are
        % saved in a vector and also every single pair is considered
        indices = nchoosek(1:obj.nExperiments,2);
        DistInput_proto =  obj.InputData(indices(:,1),:)-obj.InputData(indices(:,2),:);
        
        % Euclidean Distance
        obj.DistInput = sqrt(sum(DistInput_proto.^2,2));

        % Estimation of Variance is the half of the square of the 
        % differences of outputs 
        obj.VarianceEstimation =  ((obj.OutputData(indices(:,1),:)-obj.OutputData(indices(:,2),:)).^2)/2;
        obj.VarianceEstimation = obj.VarianceEstimation;

        % Check
        expectedSize=(obj.nExperiments)*(obj.nExperiments-1)/2;
        if size(obj.DistInput,1)~=expectedSize
            error('DistInput has dimension %ix%i but expected is 1x%i',size(obj.DistInput,1),size(obj.DistInput,2),expectedSize);
        end

        if size(obj.VarianceEstimation,1)~=expectedSize
            error('VarianceEstimation has dimension %ix%i but expected is 1x%i',size(obj.VarianceEstimation,1),size(obj.VarianceEstimation,2),expectedSize);
        end
    elseif obj.CovariogramUsesAbsoluteDistance
        % DistInput is not anymore just a vector but a matrix of
        % dimension (nExperiments*(nExperiments-1)/2)XnInputVar. Each row
        % represents a pair of input values which are used for the
        % OutputData. Each column contains the abolute difference of
        % one input variable between the difference measurements

        % Calculate for every combination of the experiments
        try
            indices = VChooseK(1:obj.nExperiments,2);
        catch ex
            warning('VChooseK could not be used. Use Matlabs nchoosek insteads')
            warning(ex.message);
            indices = nchoosek(1:obj.nExperiments,2);
        end
        
        obj.DistInput =  abs(obj.InputData(indices(:,1),:)-obj.InputData(indices(:,2),:));
        obj.VarianceEstimation =  ((obj.OutputData(indices(:,1),:)-obj.OutputData(indices(:,2),:)).^2)/2;
    elseif obj.CovariogramUsesAbsoluteValues
    else
        error('Either "CovariogramUsesEuclideanDistance" or "CovariogramUsesAbsoluteDistance" or "CovariogramUsesAbsoluteValue" have to be "true"')
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
