function [] = estimateBasisFctParametersFminCon(obj)
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.
    %% Calculate
    obj.covariance_zero = [];
    
    % Initialization
    if isempty(obj.InitialBasisFctParameters)
        obj.InitialBasisFctParameters = ones(obj.nBasisFctParameters,1);
    end
    
    % Actual Optimization
    if isempty(obj.OptimizationOption)
        options = obj.setOptionsFminCon();
    else
        options = obj.OptimizationOption;
    end
    
    if obj.UseAuxiliaryCondition
        [optPara, minError] ...
         = fmincon(@obj.calcKrigingParaEstimationObjFct,...
                   obj.InitialBasisFctParameters,...
                   [],[],[],[],...
                   obj.LBBasisFctParameters,obj.UBBasisFctParameters,...
                   @calcKrigingParaEstimationCondition,options);
    else
        [optPara, minError] ...
         = fmincon(@obj.calcKrigingParaEstimationObjFct,...
                   obj.InitialBasisFctParameters,...
                   [],[],[],[],...
                   obj.LBBasisFctParameters,obj.UBBasisFctParameters,...
                   [],options);
    end
    
    obj.BasisFctParameters = optPara;
    obj.ModelErrorBasisFct = minError;

    %% Nested Function
    function [saveCond,cEq] = calcKrigingParaEstimationCondition(para)
    % Calculates the quality of the parameter estimation with respect to the
    % Kriging quality measurement: inv(C)*Z=!0
        cEq =0;
        % Create extended Covariance Matrix
            % Initialization
        if isempty(obj.getCovariogramMatrix)
            obj.calcCovariogramMatrix
        end
        CovMatrix = obj.getCovariogramMatrix;
        quality = 0;

        if length(para)~=obj.getnBasisFctParameters
            warning('Number of given parameters (%i) and the needed number of parameters (%i) differes from each other',length(para),obj.getnBasisFctParameters)
        end
        saveCond = zeros(obj.nBasisFctDerivative,1);
        % Separate the different 
        C = zeros(obj.getnExperiments+1+1);
        C(1:obj.getnExperiments,1:obj.getnExperiments) = ...
                            CovMatrix(1:obj.getnExperiments,1:obj.getnExperiments);
        for iDerivative=1:obj.nBasisFctDerivative
            % Basis Function
            C(1:obj.getnExperiments,obj.getnExperiments+1) = obj.getBasisFct{1}(para,obj.getInputData);
            % Derivatives
            C(1:obj.getnExperiments,end) = obj.getBasisFctDerivative{iDerivative}(para,obj.getInputData);
            % All together
            C=triu(C)+triu(C)'-bsxfun(@times,eye(length(C)),diag(C));

            % Calculate the objective function
            Cinv = inv(C);
    %         saveCoeff(iDerivative) =(Cinv(end,:)*[obj.getOutputData;zeros(2,1)])/sum(Cinv(end,:).^2);
            saveCond(iDerivative) = sum(Cinv(end,:).^2);
        end
        saveCond = -saveCond+10;
    %     quality=sum(saveCoeff.^2);

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
