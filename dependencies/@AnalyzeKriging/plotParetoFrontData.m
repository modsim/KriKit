function [] = plotParetoFrontData(obj,varargin)
% []=plotParetoFrontData(krigingObjectIndex)
%
% This function allows to plot the data in the obective space and indicates
% Pareto Optimal sample points with a different color. 
%
% Input: 
% krigingObjectIndex ... vector which contains the indices of the
%                        objectives of interest. 
% 
% You can set: 
% - ShowParetoRectangles ... if true, Pareto optimal points are connected
%                            by rectangles, indicating pareto optimal
%                            regions. (Only 2D case) 
% - WriteNumbers ... Pareto optimal points are numerated with increasing
%                    values in the first obective.
% 
% You can get: -
%
%
% Note: Run "determineParetoSet()" before using "plotParetoFrontData()".
% This function does not open an extra figure but uses the current one
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

% Initialization
krigingObjectIndex = varargin{1};
[nObj,outputMatrix] = initializeAndCheckInput(krigingObjectIndex);


hold on
[HandlingSamples,HandlingSamplesOptimal,paretoSamples]=sortAndPlot(outputMatrix,[158,238,255;45,202,176]/255);

% Create Legend
if isempty(HandlingSamples)
    legend(HandlingSamplesOptimal,'Optimal Sample Points')
else
    legend([HandlingSamples,HandlingSamplesOptimal],'Non-optimal Sample Points','Optimal Sample Points')
end

% Numerate Samples if desired
if obj.WriteNumbers==1
    numerateSamples(nObj,paretoSamples);
end

%% Nested Function
% -------------------------------------------------------------------------
function [nObj,outputMatrix] = initializeAndCheckInput(krigingObjectIndex)
    nObj = length(krigingObjectIndex);
    if ~(length(krigingObjectIndex)==2||length(krigingObjectIndex)==3)
        error('plotParetoFrontData can only be used for two objectives but "krigingObjectIndex" has the length %i',length(krigingObjectIndex))
    end

    outputMatrix = zeros(obj.KrigingObjects{krigingObjectIndex(1)}.getnExperiments,nObj);
    for iObjNested = 1:nObj
        outputMatrix(:,iObjNested) = obj.KrigingObjects{krigingObjectIndex(iObjNested)}.getOutputData;
    end
end
% -------------------------------------------------------------------------
function []=numerateSamples(nObj,paretoSamples)
    for iParetoPoint = 1:obj.nParetoSetExperiments
        switch nObj
            case 2
                text(paretoSamples(iParetoPoint,1),...
                    paretoSamples(iParetoPoint,2),...
                    num2str(iParetoPoint),...
                    'HorizontalAlignment','center','FontWeight','bold');
            case 3
                text(paretoSamples(iParetoPoint,1),...
                    paretoSamples(iParetoPoint,2),...
                    paretoSamples(iParetoPoint,3),...
                    num2str(iParetoPoint),...
                    'HorizontalAlignment','center','FontWeight','bold');
        end
    end
end
% -------------------------------------------------------------------------
function [HandlingSamples,HandlingSamplesOptimal,paretoSamples]=sortAndPlot(data,colorVec)
    % Size of symbol indicating the sample points
    markerSize = 7.5;
    
    % If not pareto set is determined
    if obj.nParetoSetExperiments<=0
        warning('No Pareto Set was determined. "determineParetoSet" was called')
        determineParetoSet(obj,krigingObjectIndex)
    end
    
    % Sort ParetoSet allong th x-axis
    paretoSamples = sortrows(obj.ParetoSetExperiments);
    
    % Connect pareto points by a "stair" indicating dominated and non
    % domiated are
    if obj.ShowParetoRectangles~=false&&nObj==2
        [ExtSamples1,ExtSamples2]=manhattan(paretoSamples(:,1),paretoSamples(:,2),0);
        plot(ExtSamples1,ExtSamples2,'k');
        [ExtSamples1,ExtSamples2]=manhattan(paretoSamples(end:-1:1,1),paretoSamples(end:-1:1,2),0);
        plot(ExtSamples1,ExtSamples2,'k--');
    end
    
    % Plot data
    if size(paretoSamples,1)<size(data,1)
        switch nObj
            case 2
                HandlingSamples = plot(data(:,1),data(:,2),'kO','MarkerFaceColor',colorVec(1,:),'MarkerSize',markerSize);
            case 3
                HandlingSamples = plot3(data(:,1),data(:,2),data(:,3),'kO','MarkerFaceColor',colorVec(1,:),'MarkerSize',markerSize);
            otherwise
        end
    else
        HandlingSamples = [];
    end
    
    % Indicate Pareto optimal points by different color
    switch nObj
        case 2
            HandlingSamplesOptimal = plot(paretoSamples(:,1),paretoSamples(:,2),'kO','MarkerFaceColor',colorVec(2,:),'MarkerSize',markerSize);
        case 3
            HandlingSamplesOptimal = plot3(paretoSamples(:,1),paretoSamples(:,2),paretoSamples(:,3),'kO','MarkerFaceColor',colorVec(2,:),'MarkerSize',markerSize);
        otherwise
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
