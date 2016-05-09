function [predictionMatrix] = doMututalPrediction(obj,varargin)
% [predictionMatrix] = doMututalPrediction(indicesOfKrigingObjects,inputValues)
%
% You can set: -
% - ScaleFactorMutalPrediction ... adjust smoothness for the transient region
%                                between Krigign objects.
%                                0<=ScaleFactorMutalPrediction<=1, 
%                                ScaleFactorMutalPrediction=0 boolean
%                                transition
% - LBInputVarInterpolation, UBInputVarInterpolation ... Defines the
%       range for each Kriging object in which it is used for the
%       interpolation
%
% You can get: - 
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

%% Initialization
indicesOfKrigingObjects = varargin{1};
inputValues = varargin{2};
nChosenKrigingObjects = length(indicesOfKrigingObjects);
nInputVar = obj.KrigingObjects{indicesOfKrigingObjects(1)}.getnInputVar;

% Normalize dimensions of indicesOfKrigingObjects
if ~any(size(indicesOfKrigingObjects)==1)
    error('indicesOfKrigingObjects must be an array')
end
if size(indicesOfKrigingObjects,1)==1&&size(indicesOfKrigingObjects,2)>1
    indicesOfKrigingObjects = indicesOfKrigingObjects';
end
if size(inputValues,1)==nInputVar&&size(inputValues,2)~=nInputVar
    inputValues = inputValues';
end

% Check Content
indicesOfKrigingObjects = round(indicesOfKrigingObjects);
if any(indicesOfKrigingObjects<0|indicesOfKrigingObjects>obj.getnKrigingObjects)
    error('Only indices between 0 and maximal number of Kriging objects (%i) are allowed',obj.getnKrigingObjects)
end

%% Prework
% Identify bound of different Krigigng objects
LBobjects = zeros(nChosenKrigingObjects,nInputVar);
UBobjects = zeros(nChosenKrigingObjects,nInputVar);
for iKrigingObject = 1:nChosenKrigingObjects
    indexChoosen = indicesOfKrigingObjects(iKrigingObject);
    if nInputVar~=obj.KrigingObjects{indexChoosen}.getnInputVar
        error('Number of input variables of chosen Kriging objects must be the same. Number of Input variables of Object %i and Object %i are not the same',1,iKrigingObject)
    end
    defineBoundOfInputVar(obj,indicesOfKrigingObjects(iKrigingObject))
    
    LBobjects(iKrigingObject,:)= obj.LBInputVarInterpolation{indexChoosen};
    UBobjects(iKrigingObject,:)= obj.UBInputVarInterpolation{indexChoosen};
    
end

% Calculate weights Weights
weights = obj.calcWeightsForMutualPrediciton(inputValues,LBobjects,UBobjects);

%% Actual Prediction

% Inndividual Predictions
predictionMatrix = zeros(size(inputValues,1),2*nChosenKrigingObjects);
for iKrigingObject = 1:nChosenKrigingObjects
    predictionMatrix(:,(iKrigingObject-1)*2+1:2*iKrigingObject)=obj.KrigingObjects{indicesOfKrigingObjects(iKrigingObject)}.prediction(inputValues);
end

% Do Smoothing
weightedAverage=sum(predictionMatrix(:,1:2:end).*weights,2)./(sum(weights,2));
sigmaWeighted = sum(predictionMatrix(:,2:2:end).*weights.^2,2)./(sum(weights,2).^2); 

% Final Output
predictionMatrix = [weightedAverage,sigmaWeighted];
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
