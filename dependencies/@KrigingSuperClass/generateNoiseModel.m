function [] = generateNoiseModel(obj,varargin)
% [] = generateNoiseModel(sigmaError)
%
% This function generate a model for the output noise heteroscedastic. The
% model is inheritly used for the actual prediction .
%
% Input: 
% - sigmaError ... Initial measurement stadard deviation which is used for
%                  the as homosedatic noise as starting point. If not
%                  provided,sigmaError is set the value estimated by the
%                  maximum likelihood estimator
% 
% Output:
%
% You can set: 
% - ThresholdSolver ... optimization stops as soon the relative difference
%                       (maxDiffPred) is lower than ThresholdSolver
% - nIterationsSolver ... optimization stop if nIterationsSolver iterations
%                         were run
% - CovariogramModelChoice ... Covariogram which is used for the estimation
%                              of the output noise
%
% You can get:
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

% if obj.NormOutput
%     error('Not yet implemented for normalized output')
% end

% Reset Noise Model
obj.KriKitObjNoise = OrdinaryKriging;
obj.KriKitObjNoise.setMaxSizeOfPredictions(obj.getMaxSizeOfPredictions)
obj.KriKitObjNoise.setShowWaitingBar(false);
UseMatlabRegressionGPBackUp = obj.UseMatlabRegressionGP;
obj.UseMatlabRegressionGP = false;
sigmaErrorBackUp = obj.sigmaError;
if ~isempty(varargin)
    obj.sigmaError = varargin{1};
    obj.calcCovariogramMatrix
end

% Generate initial Noise set
if obj.ShowDetails
    fprintf('Generating Noise Model\n');
end
predictionData = obj.prediction(obj.getInputData);
nSamples = 1e3;
predictedSamples = randn(obj.nExperiments,nSamples);
predictedSamples = bsxfun(@times,predictedSamples,predictionData(:,2));
predictedSamples = bsxfun(@plus,predictedSamples,predictionData(:,1));
arithmeticMean = bsxfun(@minus,predictedSamples,obj.getOutputData);
varDelta = 1/(2*nSamples)*sum((arithmeticMean).^2,2);

% Run Optimization
for iIter = 1 : obj.nIterationsSolver

    % Set Up Noise Model
    obj.KriKitObjNoise.setInputData(obj.getInputData)
    obj.KriKitObjNoise.setOutputData(log(varDelta))
    obj.KriKitObjNoise.setCovariogramModelChoice(obj.CovariogramModelChoice);
    obj.KriKitObjNoise.generateRegressionGPModel
    predictionData2 = obj.prediction(obj.getInputData);
    
    % Estimate variation
    nSamples = 1e3;
    predictedSamples = randn(obj.nExperiments,1e3);
    predictedSamples = bsxfun(@times,predictedSamples,predictionData(:,2));
    predictedSamples = bsxfun(@plus,predictedSamples,predictionData(:,1));
    arithmeticMean = bsxfun(@minus,predictedSamples,obj.getOutputData);
    varDelta2 = 1/(2*nSamples)*sum((arithmeticMean).^2,2);
    
    % Check if converged
    maxDiffPred = max(abs(predictionData(:,2)-predictionData2(:,2))./predictionData(:,2));
    maxDiff = max(abs(sqrt(varDelta2(:,1))-sqrt(varDelta(:,1))));
    if obj.ShowDetails
        fprintf('iIter = %i - maxDiff = %e - maxDiffPred = %e \n',iIter,maxDiff,maxDiffPred);
    end
    if maxDiffPred<obj.ThresholdSolver
        fprintf('Optimization has converged\n');
        break
    end
    
    % Update
    predictionData = predictionData2;
    varDelta = varDelta2;
end
if iIter==obj.nIterationsSolver
    fprintf('Optimization stopped as maximimum iteration steps are reached \n');
end
predVariance = (obj.KriKitObjNoise.prediction(obj.getInputData));

% Save Results and recalculate the covariogram matrix
obj.setHeterogeneousNoise(sqrt(exp(predVariance(:,1))))
obj.calcCovariogramMatrix

% Restore Original values
obj.UseMatlabRegressionGP = UseMatlabRegressionGPBackUp;
if ~isempty(varargin)
    obj.sigmaError = sigmaErrorBackUp;
    obj.calcCovariogramMatrix
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
