function [weights] = calcWeightsForMutualPrediciton(obj,pointOfInterest,lb,ub)
% [weights] = calcWeightsForMutualPrediciton(pointOfInterest,lb,ub)
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

% Calculates the weights used for the weighted sum of Kriging predictions
% in doMutualPrediction(). If the point of interest is in the defined range
% for the kriging objective, the weight value is set to 1, otherwise weight
% value is descreasing with the distance to the defined range:
%           weight = min(1,1/(distanceToRange)*ScaleFactorMutalPrediction)
%
%
% Input:
% pointOfInterest ... matrix containting the interpolation points
%                     (nInterpolationPointsXnInputVar) 
% lb, ub ... define the border of the Interpolation ranges
%            (nKrigingObjectsXnInputVar)
%
% You can set:
% ScaleFactorMutalPrediction ... adjust smoothness when point of interest
%                                is not in the defined range
%
% You can get: - 
%
% Copyright 2014-2015: Lars Freier, Eric von Lieres
% See the license note at the end of the file.
    

if sum(size(lb)==size(ub))<2
    error('Dimensions of lb and ub must be the same')
end

% scaleFactor = max(min(min(ub-lb))/5e4);
scaleFactor = obj.ScaleFactorMutalPrediction;

weights = zeros(size(pointOfInterest,1),size(lb,1));
for iPoint = 1:size(pointOfInterest,1)
    differenceLB = bsxfun(@minus,pointOfInterest(iPoint,:),lb);
    differenceUB = bsxfun(@minus,pointOfInterest(iPoint,:),ub);
    for iSection = 1:size(lb,1)
        
        % Points inside the lb and ub defined range get the weight 1.
        % Otherwise, weight is discreasing with increasing distance to
        % defined range
        if sum(differenceLB(iSection,:)>=0)~=size(lb,2)||sum(differenceUB(iSection,:)<=0)~=size(ub,2)
            % difference to the lower bound (sum twice since otherwise error if diff is empty)
            a = sum(sum(differenceLB(iSection,differenceLB(iSection,:)<0).^2));
            b = sum(sum(differenceUB(iSection,differenceUB(iSection,:)>0).^2));
            weights(iPoint,iSection) = min(1/((a+b)/scaleFactor),1);
        else
            weights(iPoint,iSection) = 1;
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
