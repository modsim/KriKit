function [] = findPotentialOutlier(obj,varargin)
% [] = findPotentialOutlier(obj,KrigingIndex)
%
% Determines Outliers. As outliers lie outside of the confidence interval.
% The confidence interval is defined by "WidthConfidenceInterval"
%
% You can set:
% - WidthConfidenceInterval ... defined width of the confidence interval
%                               (mean +- WidthConfidenceInterval*sd) 
%
% You can get: - 
% PotentialOutlier ... cell where PotentialOutlier{KrigingIndex} is a
%                      matrix containing the ouliers (nOutliersXnInputVar) 
% PotentialToBeOutlier ... cell where PotentialToBeOutlier{KrigingIndex} is a
%                      matrix containing the relative distance to the
%                      Kriging prediction(=distacne/predictionValue)
%                      (nOutliersX1) 
% PotentialOutlierBinary ... vector which entries are true if the sample
%                            point is an outlier (nOutliersX1)
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.
    
%% Initialization
KrigingIndex = varargin{1};
if (size(KrigingIndex)==1)~=1
    error('KrigingIndex must be a scalar')
end

%% Prediction
inputData = obj.KrigingObjects{KrigingIndex}.getInputData;
prediction = obj.KrigingObjects{KrigingIndex}.prediction(inputData);

%% Check if outside of confidence interval
diffPredvsOutput = abs(prediction(:,1)-obj.KrigingObjects{KrigingIndex}.getOutputData);
obj.PotentialOutlierBinary=diffPredvsOutput>=obj.WidthConfidenceInterval*prediction(:,2);

% Normalize and sort 
RelativeDiffPredvsOutput = diffPredvsOutput(obj.PotentialOutlierBinary)./prediction(obj.PotentialOutlierBinary,2);
inputDataOutlier = inputData(obj.PotentialOutlierBinary,:);
[diffPredvsOutput_sort,indicesData] = sort(RelativeDiffPredvsOutput,'descend');

% Save
obj.PotentialOutlier{KrigingIndex} = inputDataOutlier(indicesData,:);
obj.PotentialToBeOutlier{KrigingIndex} = diffPredvsOutput_sort;

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
