function [ExpectedImprovement]=calcExpectedImprovement_MO(obj,varargin)
% [ExpectedImprovement]=calcExpectedImprovement_MO(objectiveIndices,estimationPoints)
% 
% This function calculates the exprective improvement for the
% multi-objective case. Improvement to an increasing in the
% hypervolume(EIHV). Basic concepts of EIHV are described in Emmerich et
% al. 2006 "Single- and Multiobjective Evolutionary Optimization Assisted
% by Gaussian Random Field Metamodels"
% 
% Input:
% - objectiveIndices ... indices of kriging objects of interest [1XnObj]
% - estimationPoints ... Matrix of input points [nPointsXnInputVar]
%
% Output:
% - ExpectedImprovement ... veector containing the expected improvment value
%                           associated with input points given in
%                           "estimationPoints" [nPointsX1]
%
% You can set:
% - InequalityConstraintOutputHandle ... (only in case of BayesianOptimization): 
%           A function handle of the form
%           [boolConstraintHolds]=inequalityConstraintOutput(estimationMean,estimationSD)
%           With 
%           estimationMean ... kriging prediction at the point of interest
%                              given in "estimationPoints" [nPointsX1]
%           estimationSD ... kriging prediction error at the point of
%                            interest, given in "estimationPoints"
%                            [nPointsX1] 
%           boolConstraintHolds ... boolean output vector containing true
%                                   if the point hold the constraint
%                                   [nPointsX1] 
%
% - DegreeOfExpectedImprovement ... exploration factor for adjusting the
%                   trade-off between exploitation and exploration. 
%                   DegreeOfExpectedImprovement=0: Only Kriging prediction
%                       is considered --> Local optimization
%                   DegreeOfExpectedImprovement->inf: Only uncertain
%                       regions are considered --> gloabl optimization
% - ReferencePointHyperVolume ... Reference point used for calculating the
%                                 hypervolume
%
% You can get: 
%
% Note: Before Running "calcExpectedImprovement_MO" you should calculate
% the run "determineParetoSet()"
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

%% Initialization
% Identify Input
objectiveIndices = varargin{1};
estimationPoints = varargin{2};

if size(estimationPoints,2)~=obj.KrigingObjects{objectiveIndices(1)}.getnInputVar
    error('estimationPoints must have nInputVar(=%i) columns',obj.KrigingingObjects{objectiveIndices(1)}.getnInputVar)
end

% Initialization
nObj = length(objectiveIndices);
nEstimationPoint = size(estimationPoints,1);
estimationMean = zeros(nEstimationPoint,nObj);
estimationSD = zeros(nEstimationPoint,nObj);

%% Do prediction and adjust exploration behaviour by modifying prediction error
DegreeOfExpectedImprovement = obj.getDegreeOfExpectedImprovement;
for iObj=1:nObj
    krigingPrediction = obj.KrigingObjects{objectiveIndices(iObj)}.prediction(estimationPoints);
    estimationMean(:,iObj) = krigingPrediction(:,1);
    estimationSD(:,iObj) = DegreeOfExpectedImprovement*krigingPrediction(:,2);
end

% Check if inequality constraits for the output hold
if isprop(obj,'InequalityConstraintOutputHandle')&&~isempty(obj.getInequalityConstraintOutputHandle)
    % InequalityConstraintOutputHandle is part of BayseanOptimization not
    % AnalyzeKriging
    handleInEq = obj.getInequalityConstraintOutputHandle;
    boolIsValid = handleInEq(estimationMean,estimationSD);
else
    boolIsValid = true(nEstimationPoint,1);
end

%% Determine Pareto set
% Convert problem to a minimization problem
if isempty(obj.getParetoSetExperiments)||size(obj.getParetoSetExperiments,2)~=nObj
    warning('ParetoSetExperiments has to be determined before calcExpectedImprovement_MO can be used')
    obj.determineParetoSet(objectiveIndices)
end
paretoSet = bsxfun(@times,obj.getParetoSetExperiments,-obj.MinMax(objectiveIndices));
estimationMean = bsxfun(@times,estimationMean,-obj.MinMax(objectiveIndices));

if isempty(obj.ReferencePointHyperVolume)
    refPoint = max(paretoSet);
else
    refPoint = obj.ReferencePointHyperVolume;
end
refPoint = bsxfun(@times,refPoint,-obj.MinMax(objectiveIndices));
%% Actual Calculation
ExpectedImprovement = zeros(nEstimationPoint,1);
switch nObj
    case 2
        for iEstimationPoint = 1:size(estimationMean,1)
            if boolIsValid(iEstimationPoint)
                % Implementation based on the work of Emmerich et al. 2011
                % "Hypervolume-based Expected Improvement: Monotonicity
                % Properties and Exact Computation"
                ExpectedImprovement(iEstimationPoint) = calcEHIV_2D(paretoSet,refPoint,...
                                                                    estimationMean(iEstimationPoint,:),...
                                                                    estimationSD(iEstimationPoint,:));
            else
                ExpectedImprovement(iEstimationPoint) = 1e-10;
            end
        end
    case 3
        for iEstimationPoint = 1:size(estimationMean,1)
            if boolIsValid(iEstimationPoint)
                % Implementation is based on work of Hubkens et al. 2015
                % "Faster Exact Algorithms for Computing Expected
                % Hypervolume Improvement"
                ExpectedImprovement(iEstimationPoint) = calcEIHVmex(-paretoSet,...
                                                                   -refPoint,[-estimationMean(iEstimationPoint,:),...
                                                                   estimationSD(iEstimationPoint,:)]);
            else
                ExpectedImprovement(iEstimationPoint) = 1e-10;
            end
        end
    otherwise
        for iEstimationPoint = 1:size(estimationMean,1)
            if boolIsValid(iEstimationPoint)
                % Expected improvement is calculated via well defined
                % markov chain monte carlo chain approach
                ExpectedImprovement(iEstimationPoint) = calcEIMO_MCMC(obj,-paretoSet,...
                                                                          -refPoint,-estimationMean(iEstimationPoint,:),...
                                                                          estimationSD(iEstimationPoint,:));
            else
                ExpectedImprovement(iEstimationPoint) = 1e-10;
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
