function [PearsonCoeff] = plotOutputvsPrediction(obj,varargin)
% [PearsonCoeff] = plotOutputvsPrediction(KrigingObjectIndex)
%
% This fucntion creates a 2D-plot which allows to evaluate the Kriging
% model by plotting the Kriging prediction with the actual output.
% Approximatly, 67% of all data point should be inside the errorbar.
%
% Input: 
% - KrigingObjectIndex ... index of kriging objects of interest [1X1]
%
% Output:
% - PearsonCoeff ...Perason correlation coefficient between output and
%                   Kriging prediction
%
% You can set:
% - ShowBounds ... if true error bars are used. If false, output data are
%                  plotted agains Kriging prediction
%
% You can get:
% - OutsideOfConfi ...
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

KrigingObjectIndex = varargin{1};
if length(KrigingObjectIndex)~=1
    error('KrigingObjectIndex has to be a scalar')
end

% Do Prediction
predictionOutput = obj.KrigingObjects{KrigingObjectIndex}.prediction(...
                        obj.KrigingObjects{KrigingObjectIndex}.getInputData);

% Make Plot 
% 'ko','MarkerFaceColor',[255/255 102/255 0/255]);
figure();
hold on;
if obj.ShowBounds==1
    errorbar(obj.KrigingObjects{KrigingObjectIndex}.getOutputData,predictionOutput(:,1),...
        predictionOutput(:,2),'ko','MarkerFaceColor',[255/255 102/255 0/255]);
else
    plot(obj.KrigingObjects{KrigingObjectIndex}.getOutputData,predictionOutput(:,1),...
        'ko','MarkerFaceColor',[255/255 102/255 0/255]);
end
plot(obj.KrigingObjects{KrigingObjectIndex}.getOutputData,...
    obj.KrigingObjects{KrigingObjectIndex}.getOutputData,'k')
legend('Output vs Prediction','Reference Line','Location','NorthWest')
xlabel('Output','FontSize',20)
ylabel('Prediction','FontSize',20)
title(strcat('Output vs Prediction (',obj.KrigingObjectNames{KrigingObjectIndex},')'),'FontSize',20)

%% Pearson Coefficient
pearsonCoeffMatrix = corrcoef(obj.KrigingObjects{KrigingObjectIndex}.getOutputData,predictionOutput(:,1));
PearsonCoeff = pearsonCoeffMatrix(1,2);

%% Relative Amount of Points outside the confidence interval
UB = predictionOutput(:,1) + predictionOutput(:,2);
LB = predictionOutput(:,1) - predictionOutput(:,2);
% Points Outside of the Upper bound of the confi-interval
biggerThanUpperBound = sum((obj.KrigingObjects{KrigingObjectIndex}.getOutputData-UB)>0);
lowerThanLowerBound = sum((obj.KrigingObjects{KrigingObjectIndex}.getOutputData-LB)<0);

obj.OutsideOfConfi = (biggerThanUpperBound + lowerThanLowerBound)/...
                            obj.KrigingObjects{KrigingObjectIndex}.getnExperiments;
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
