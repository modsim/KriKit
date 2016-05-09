function [ExpectedImprovement] = calcExpectedImprovementMainPart(obj,varargin)
% [ExpectedImprovement] = calcExpectedImprovementMainPart(KrigingObjectIndex,prediction)
%
% Input:
% - KrigingObjectIndex ...  Index of Kriging object of interest
% - prediction ... prediction made by Kriging at the point of interest
%
% Output:
% - ExpectedImprovement ... expected improvement at the point of interest
%                           [nPointsX1]
%
%
% You can set: -
%
% You can get: -
%
% For further optional settings see "calcFinalExpImprovement()"
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

%% Initialization
KrigingObjectIndex = varargin{1};
prediction = varargin{2};

if length(KrigingObjectIndex)>1
    error('calcExpectedImprovement can not be used for multi-objective or mutual Kriging')
end

% Convert Output to a minimization problem
if length(varargin)>2
    checkPoints = varargin{3};
    bestOutput = -obj.MinMax(KrigingObjectIndex)*obj.KrigingObjects{KrigingObjectIndex}.prediction(checkPoints);
    bestOutput = bestOutput(:,1);
else
    % Calculate the expected improvements at the particular sample
    % locations
    outputData = -obj.MinMax(KrigingObjectIndex)*obj.KrigingObjects{KrigingObjectIndex}.getOutputData;
    [~,indexMin] = min(outputData);
    bestOutput = outputData(indexMin);
end


% Argument for cumulative distribution function
u = (bestOutput - prediction(:,1))./prediction(:,2);

%% Calculate expected improvement
[ExpectedImprovement] = obj.calcFinalExpImprovement(prediction,u);

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
