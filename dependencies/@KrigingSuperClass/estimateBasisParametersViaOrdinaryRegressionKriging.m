function [] = estimateBasisParametersViaOrdinaryRegressionKriging(obj)
% [] = estimateBasisParametersViaOrdinaryRegressionKriging()
%
% This function does regression kriging in order to estimate the nonlinear
% parameter of the basis function. Here, not the measurements itself are
% used for the LSQ regression but, the approximation using ordinary Kriging.
% This procedure might be more robust to outliers.
%
% Note: internally a new kriging model is calculated. The current one is
% not changed
%
%
% You can set:
% - BasisFct
% - CovariogramModelChoice
%
% You can get:
% - OrdinaryRegressionKrigingObject
% - BasisFctParameters
% - ModelErrorVar
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

% ##### Make checks #####
if obj.nBasisFct~=1
    error('estimateBasisFctParametersViaOrdinaryRegressionKriging is only defined for one basis function but you use %i basis functions',obj.nBasisFct)
end
if obj.nBasisFctParameters<1
    error('Basis function does not have any parameters')
end

% ##### Use ordinary kriging to calculate the interpolation near the measurement points #####
doOrdinaryKrigingInterpolation;

% ##### Do non-linear regression #####
doRegressionGA()


%% ------------------------------------------------------------------------
    function [] = doOrdinaryKrigingInterpolation()
            % Declare Kriging object
        obj.OrdinaryRegressionKrigingObject = OrdinaryKriging;
        obj.OrdinaryRegressionKrigingObject.setCovariogramEstimationType(obj.getCovariogramEstimationType())
        obj.OrdinaryRegressionKrigingObject.setNormInput(obj.getNormInput);
        obj.OrdinaryRegressionKrigingObject.setNormOutput(obj.getNormOutput);
            % Copy solver options
        obj.OrdinaryRegressionKrigingObject.setUseSolver(obj.getUseSolver)
        obj.OrdinaryRegressionKrigingObject.setShowDetails(obj.getShowDetails)
            % Load data
        obj.OrdinaryRegressionKrigingObject.setInputData(obj.getInputData);
        obj.OrdinaryRegressionKrigingObject.setOutputData(obj.getOutputData);
            % Choose Covariogram model
        if obj.CovariogramModelChoice<0
            obj.OrdinaryRegressionKrigingObject.setCovariogramModelChoice(2);
        else
            obj.OrdinaryRegressionKrigingObject.setCovariogramModelChoice(...
                obj.getCovariogramModelChoice());
        end
        obj.OrdinaryRegressionKrigingObject.setUseMatlabRegressionGP(obj.getUseMatlabRegressionGP);
%             % estimate variogram
        estimateVariogram
            % Estimate the output variables near the measuremnt points. Not
            % excactly at the measurment points since here the nugged has no effect
        obj.EstimatedOutputDataViaOrdinaryKriging = obj.OrdinaryRegressionKrigingObject.prediction(obj.InputData+1e-8);
        % -----------------------------------------------------------------
        function [] = estimateVariogram()
            % ##### Display ##### 
            if obj.OrdinaryRegressionKrigingObject.getShowDetails
                fprintf('Variogram estimation:\n')
            end
            % ##### Actual estimation ##### 
            obj.OrdinaryRegressionKrigingObject.makePrework;
            % ##### Show result ##### 
            if obj.OrdinaryRegressionKrigingObject.getShowDetails
                obj.OrdinaryRegressionKrigingObject.plotVariogram;
            end
        end
    end
%% ------------------------------------------------------------------------
    function [] = doRegressionGA()
    % Uses genetic algorithm to estimate the 
            % ##### Define Options #####
            switch obj.ShowDetails
                case 0
                    showDetails='off';
                case 1
                    showDetails='iter';
                    fprintf('Basis parameter estimation estimation:\n')
                otherwise
                    error('ShowDetails=%i was not defined',ShowDetails);
            end
            if obj.PopulationSize<=0
                obj.PopulationSize=min(max(10*obj.nInputVar,40),100);
            end
            optionsGA = gaoptimset('Generations',obj.Generations,'TimeLimit',120,'display',showDetails,'PopulationSize',obj.PopulationSize);
            
            % ##### DO regression #####
            [optPara, minError] = ga(@calcLeastSquareError,obj.nBasisFctParameters,...
                           [],[],[],[],...
                           [],[],...
                           [],optionsGA);

            % ##### Final output #####
            obj.BasisFctParameters = optPara;
            obj.ModelErrorVar = minError;
%             keyboard
            
        % ----------------------- Nested Function -----------------------
        function [SQ] = calcLeastSquareError(varargin)
            parameters = varargin{1};
            basisFunctionValues = obj.BasisFct{1}(parameters,obj.InputData);
            SQ = sum((basisFunctionValues-obj.EstimatedOutputDataViaOrdinaryKriging(:,1)).^2);
        end
    end
end

% figure()
% hold on
% plot(obj.getInputData,obj.getOutputData,'.');
% plot(obj.getInputData,obj.EstimatedOutputDataViaOrdinaryKriging(:,1));
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
