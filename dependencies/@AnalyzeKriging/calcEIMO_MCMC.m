function [EI] = calcEIMO_MCMC(obj,paretoSet,refPoint,estimationMean,estimationSD)
% [EI] = calcEIMO_MCMC(obj,paretoSet,refPoint,estimationMean,estimationSD)
%
% This function uses a Monte Carlo approach in order to approximate the
% expected hypervolume improvement.
%
% Input: -
%
% - paretoSet ... Matrix containing the Pareto optimal objective values.
%                 (nParetoPointsXnKrigingObjective)
% - refPoint ... Reference point used for calculating the hypervolume
% - estimationMean ... Kriging prediction at the point of interest
%                      given in "estimationPoints" [1XnKrigingObjective]
% - estimationSD ... Kriging prediction error at the point of interest,
%                    given in "estimationPoints" [1XnKrigingObjective] 
%
% Output: -
%
% You can set:
% - nMCMCLinksEI ...number of samples used for Monte Calor Algorithm
%
%
% You can get: -
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

%% Initialization
if size(estimationMean,1)>1||size(estimationSD,1)>1
    error('calcEIMO_MCMC can only considere one point at the time. \n "estimationMean" and "estimationSD" should be row vector')
end
nObj = size(estimationMean,2);
improvement = zeros(obj.nMCMCLinksEI,1);
oldHyperVol = Hypervolume_MEX(-paretoSet,-refPoint); % Convertion to a minimization problem

% Random Sampling with scaling
yN = bsxfun(@plus,bsxfun(@times,randn(obj.getnMCMCLinksEI,nObj),estimationSD),estimationMean);

% Calculate associated improvement for each random sample
extendedParetoSet = [paretoSet;zeros(1,size(paretoSet,2))];
if obj.nMCMCLinksEI<1
    error('nMCMCLinksEI has to be bigger equal to 1')
end
for iRandPoint = 1:obj.nMCMCLinksEI
    extendedParetoSet(end,:) = yN(iRandPoint,:);
    newHyperVol = Hypervolume_MEX(-extendedParetoSet,-refPoint);
    improvement(iRandPoint) = newHyperVol-oldHyperVol;
    if  improvement(iRandPoint)<-1e-10
        warning('Negative Hypervolume improvement')
    end
end

% Final calculation
EI = 1/obj.nMCMCLinksEI*sum(improvement);
    
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
