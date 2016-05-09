function [] = estimateBasisFctCoefficients(obj,varargin)
% estimateBasisFctCoefficients(estimateTaylor)
% This function estimates the coefficients c1:cn of the basis function save
% in "BasisFct": c1*BasisFct{1} + ... + cn*BasisFct{nBasisFct}. The
% coefficients are estimated using the pre-estimated covariogramMaxtrix.
% 
% NOTE: Call this function only after covariogram and non-linear
% basisfunction parameters were already calculated
%
% Input:
% - estimateTaylor ...If estimateTaylor=1 then the coefficients of the taylor
%                   expansion are calculated (only possible with Taylor
%                   Kriging) 
%
% You can set:
% - checkVariogram ... if false, avoid calculation of covariogram matrix.
%
%
% You can get: 
% - BasisFctCoefficients ... estimated basis function coeffcients
%
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.


%% Initialization
% Calculate the Covariogram Matrix if it doesn't exist already 
if (obj.checkVariogram==false)
    obj.calcCovariogramMatrix();
    obj.checkVariogram = true;
end

% Get information if the coeffients of the basis functions or of the taylor
% expansion should be estimated
if length(varargin)==1
    Estimation = varargin{1};
else
    Estimation = 0;
end

% ###### Do the estimation ###### 
switch Estimation
    case 1
        estimateCoeffientsOfTaylorExpansion
    case 0
        estimateCoeffientsOfBasisFct
    otherwise
        error('Input is not defined')
end
%% ------------------- Nested Function ---------------------------------
function [] = estimateCoeffientsOfBasisFct()
    % Calculated Kriging covariance matrix
    if (obj.checkVariogram==false)
        obj.calcCovariogramMatrix();
        obj.checkVariogram = true;
    end
    
    % Number of coefficients is equal to the number of basis
    % functions (nBasisFct)
    coefficientIndices = 1:obj.nBasisFct;
    obj.BasisFctCoefficients = zeros(obj.nBasisFct,1);
    
    basisFctEval = obj.CovariogramMatrix(1:obj.getnExperiments,obj.getnExperiments+1:end);
    Cov = obj.CovariogramMatrix(1:obj.getnExperiments,1:obj.getnExperiments);
    obj.BasisFctCoefficients = (basisFctEval'*(Cov\basisFctEval))\(basisFctEval'*(Cov\obj.getOutputData));
    
% % %     % Alternative
% % %     for iC = coefficientIndices
% % %         % Define the covariogram vectors (Right hand-site of the equation
% % %         % system). 
% % %         covargramVector = zeros(size(obj.CovariogramMatrix,1),1);
% % %         covargramVector(obj.getnExperiments+iC)=1;
% % %         weights = obj.CovariogramMatrix\covargramVector;
% % %         obj.BasisFctCoefficients(iC) = weights(1:obj.getnExperiments)'*obj.getOutputData;
% % %     end
    
    obj.checkVariogram = 0;
end
% -------------------------------------------------------------------------
function [] = estimateCoeffientsOfTaylorExpansion()
    % Estimate the coefficients of the taylor expansions
    
    % Initialization 
        % Number of coeffients which should be calculated (basis fct + its derivatives)
    coefficientIndices = 1:1+obj.nBasisFctParameters;
    switch obj.getTaylorExpansionOrder
        case 1
            coefficientIndices = 1:1+obj.nBasisFctParameters;
        case 2
            coefficientIndices = 1:1+obj.nBasisFctParameters + nchoosek(obj.nBasisFctParameters+2-1,2);
        otherwise
            error('TaylorExpansionOrder must be 1 or 2')
    end
    
    % Calculate extended covariance matrix
    ExtendedCovarianceMatrix=calcExtendedCovarianceMatrix();

    % ###### Calculate Coefficients ######
    estimateCoefficients;

        % Number 1: The extended covariance matrix is essential for the
        % calculation of the coefficients of the taylor expantion which
        % is done in "estimateCoefficients"
    function [ExtendedCovarianceMatrix] = calcExtendedCovarianceMatrix()
        
        % Initialize Covariogram
        ExtendedCovarianceMatrix = zeros(obj.nExperiments+length(coefficientIndices),obj.nExperiments+length(coefficientIndices));
        
        % Copy persistent part
        ExtendedCovarianceMatrix(1:obj.nExperiments,1:obj.nExperiments) = obj.CovariogramMatrix(1:obj.nExperiments,1:obj.nExperiments);
        ExtendedCovarianceMatrix = triu(ExtendedCovarianceMatrix); 
        
        % Extend Matrix by actual Basis Functions
        inputDataOriginal = obj.getInputData;
        basis = obj.getBasisFct{1}(obj.BasisFctParameters,inputDataOriginal(:,:));
        ExtendedCovarianceMatrix(1:obj.nExperiments,obj.nExperiments+1) = basis;
        
        % Extend Matrix by Derivative of Basis Functions
        for iP = 1:obj.nBasisFctParameters
            derivativeBasis = obj.getBasisFctDerivative{iP}(obj.BasisFctParameters,inputDataOriginal(:,:));
            ExtendedCovarianceMatrix(1:obj.nExperiments,obj.nExperiments+1+iP) = derivativeBasis;
        end
        
        % If the second order taylor expansion is needed
        if obj.getTaylorExpansionOrder==2
            % In the case of a second order taylor series expansion
            % the first added column represent the quadratic effect
            % of a single parameter (dh/(dp)^2). Afterwards, the quadratic
            % effect of the interaction of two parameters are added
            % (dh/(dp_i*dp_j))
            
            % Quadratic
            for iP = 1:obj.nBasisFctParameters
                derivativeBasis = obj.getBasisFctDerivative{(2-1)*obj.nBasisFctParameters+iP}(obj.BasisFctParameters,obj.InputData(:,:));
                ExtendedCovarianceMatrix(1:obj.nExperiments,obj.nExperiments+1+obj.nBasisFctParameters+iP) = derivativeBasis;
            end
            
            % Pairwise Derivative
            if obj.nBasisFctParameters>1
                combinations = nchoosek(1:obj.nBasisFctParameters,2);
                for iP = 1:size(combinations,1)
                    derivativeBasis = obj.getBasisFctDerivative{(2-1)*obj.nBasisFctParameters+obj.nBasisFctParameters+iP}(obj.BasisFctParameters,obj.InputData(:,:));
                    ExtendedCovarianceMatrix(1:obj.nExperiments,obj.nExperiments+1+2*obj.nBasisFctParameters+iP) = derivativeBasis;
                end
            end
        end
        
        % Final ExtendedCovarianceMatrix add upper lower triangle
        % matrix and subtract its diagonal since this identical in both
        % matrixes
        ExtendedCovarianceMatrix = ExtendedCovarianceMatrix + ...
                                   ExtendedCovarianceMatrix' - ...
                                   bsxfun(@times,eye(length(ExtendedCovarianceMatrix)),diag(ExtendedCovarianceMatrix));
    end
    % ---------------------------------------------------------------------
    function [] = estimateCoefficients()
        % Number 2: Actual calculation of the coefficients of the
        % taylor expansion
        
        % Initialization
        obj.BasisFctCoefficients = zeros(length(coefficientIndices),1);
        
        % Calculation
        for iC = coefficientIndices
                % Define the covargramVectors (Right hand-site of the equation system)
            covargramVector = zeros(size(ExtendedCovarianceMatrix,1),1);
            covargramVector(obj.getnExperiments+iC)=1;
                % Solve System
            weights = ExtendedCovarianceMatrix\covargramVector;
            obj.BasisFctCoefficients(iC) = weights(1:obj.getnExperiments)'*obj.getOutputData;
        end
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
