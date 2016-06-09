function [] = plotParetoFront(obj,varargin)
% [] = plotParetoFront(krigingObjVector,expectedParetoCurve,deviationParetoCurve,pCover,gridOriginal);
% Input: 
% - krigingObjVector(1XnObjectives): contains indices of kriging objects
% - expectedParetoCurve(nParetoPointsXnObjectives): matrix contains
%   coordinates of the expected pareto curve
% - deviationParetoCurve(nGridPointsX1): vector contains value
%   probability density function for each gridpoint, saved in gridOriginal
%   (see below)
% - pCover(nGridPointsX1): vector contains value cumulative probability
%   density function for each gridpoint, saved in gridOriginal
%   (see below)
% - gridOriginal(nGridPointsXnObjectives): coordinates of the objectives 
%
% Note: expectedParetoCurve,deviationParetoCurve,pCover, and gridOriginal
%       can be calculated using predictParetoCurve
%
% You can set: 
% - ShowBounds: Decide if Confidence bounce should be shown
% - ShowData: Decide if the data used for the kriging model shall be
%             plotted. For more information see documentation of
%             "plotParetoFrontData"
% - WidthConfidenceInterval: Decide how broad the condicence tube should
%                            be. This is related to normal
%                            distribution(mean value +-
%                            WidthConfidenceInterval*standard deviation).
%                            E.g.: a WidthConfidenceInterval=1 means that
%                            with a probability of 68% a pareto optimal
%                            point will be in the visualized confidence
%                            tube(/volume)
%
%
% For plotting the 3D Pareto surface, "gridfit"(Copyright (c) 2006, John D'Errico) tool which can be
% downloaded at
% http://www.mathworks.com/matlabcentral/fileexchange/8998-surface-fitting-using-gridfit'
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

    krigingObjVector = varargin{1};
    expectedParetoCurve = varargin{2};
    deviationParetoCurve = varargin{3};
    pCover = varargin{4};
    gridOriginal = varargin{5};
    
    nObj = length(krigingObjVector);
    nGridPointsOutput = round(size(gridOriginal,1)^(1/nObj));
    
    % Determine Unique Grid
    uniqueGridMatrix = calcUniqueGridMatrix();

    % Find Quantiles in order to define the confidence tube
    [gridToLow,gridToBig]=determinConfidenceBounds;
    
    switch nObj
        case 2
            figure()
            hold on
            contourf(uniqueGridMatrix(:,1),uniqueGridMatrix(:,2),reshape(deviationParetoCurve,nGridPointsOutput,nGridPointsOutput)');
            hExpected = plot(expectedParetoCurve(:,1),expectedParetoCurve(:,2),'k.','MarkerSize',15);
            
            if obj.getShowBounds
                % sort in ascend w.r.t to the first objective but descend
                % order w.r.t to the second objective
                gridSortedLB = sortrows(gridToLow,[1,-2].*obj.MinMax(krigingObjVector));
                gridSortedUB = sortrows(gridToBig,[1,-2].*obj.MinMax(krigingObjVector));
                
                hConfidence = plot(gridSortedLB(:,1),gridSortedLB(:,2),'w-.','MarkerSize',15);
                plot(gridSortedUB(:,1),gridSortedUB(:,2),'w-.','MarkerSize',15)
            else
                hConfidence = [];
            end
            
            if obj.getShowData
                plotParetoFrontData(obj,krigingObjVector);
                
                objectNonOpt = findobj(gca,'DisplayName','Non-optimal Sample Points');
                objectOpt = findobj(gca,'DisplayName','Optimal Sample Points');
                objectNonOptExtra = findobj(gca,'DisplayName','Non-optimal Extra Points');
                objectOptExtra = findobj(gca,'DisplayName','OptimalExtra Points');
            else
                objectOpt = [];
                objectNonOpt = [];
                objectOptExtra = [];
                objectNonOptExtra = [];
            end
        case 3
            % Check if "gridFit" is available
            if exist('gridfit','file')~=2
                error('This function need "gridfit". "gridfit" can be downloaded at http://www.mathworks.com/matlabcentral/fileexchange/8998-surface-fitting-using-gridfit')
            end
            
            objectiveMatrix = zeros(obj.Accuracy,nObj);
            maxGrid = max(gridOriginal);
            minGrid = min(gridOriginal);
            
            for iObj=1:nObj
                krigingIndex = krigingObjVector(iObj);
                defineBoundOfInputVar(obj,krigingIndex);
                objectiveMatrix(:,iObj) = linspace(minGrid(iObj),maxGrid(iObj),obj.Accuracy)';
            end
            
            gridObjectiveExpectedGrid = gridfit(expectedParetoCurve(:,1),expectedParetoCurve(:,2),expectedParetoCurve(:,3),objectiveMatrix(:,1),objectiveMatrix(:,2));
            
            figure
            hold on
            hExpected=surf(objectiveMatrix(:,1),objectiveMatrix(:,2),gridObjectiveExpectedGrid);
            grid on
            alpha(obj.AlphaValue)
            
            if obj.getShowBounds
                % sort in ascend w.r.t to the first objective but descend
                % order w.r.t to the remaining objective
                gridSortedLB = sortrows(gridToLow,[1,-2,-3].*obj.MinMax(krigingObjVector));
                gridSortedUB = sortrows(gridToBig,[1,-2,-3].*obj.MinMax(krigingObjVector));
                
                gridObjectiveLB = gridfit(gridSortedLB(:,1),gridSortedLB(:,2),gridSortedLB(:,3),objectiveMatrix(:,1), objectiveMatrix(:,2));
                gridObjectiveUB = gridfit(gridSortedUB(:,1),gridSortedUB(:,2),gridSortedUB(:,3),objectiveMatrix(:,1), objectiveMatrix(:,2));
                
                surf(objectiveMatrix(:,1),objectiveMatrix(:,2),gridObjectiveLB);
                hConfidence = [];
                surf(objectiveMatrix(:,1),objectiveMatrix(:,2),gridObjectiveUB);
            else
                hConfidence = [];
            end
            
            if obj.getShowData
                plotParetoFrontData(obj,krigingObjVector);
                
                objectNonOpt = findobj(gca,'DisplayName','Non-optimal Sample Points');
                objectOpt = findobj(gca,'DisplayName','Optimal Sample Points');
                objectNonOptExtra = findobj(gca,'DisplayName','Non-optimal Extra Points');
                objectOptExtra = findobj(gca,'DisplayName','OptimalExtra Points');
            else
                objectOpt = [];
                objectNonOpt = [];
                objectOptExtra = [];
                objectNonOptExtra = [];
            end
        otherwise
            error('plotParetoFront is not implemented for more than 3 objectives')
    end
    
    %% Final Adjustments
    addLegend()
    labelAxis(krigingObjVector);
    % Set Colormap
    colormap(gcf,obj.ColormapToolbox);

%% Nested
    function [uniqueGridMatrix] = calcUniqueGridMatrix()
        uniqueGridMatrix = zeros(nGridPointsOutput,nObj);
        % Bring Data Set in correct order
        for iObjNested=1:nObj
            krigingIndexNested = krigingObjVector(iObjNested);
            if obj.getMinMax(krigingIndexNested)==1
                uniqueGridMatrix(:,iObjNested)=sort(unique(gridOriginal(1:end,krigingIndexNested)),'ascend');
            else
                uniqueGridMatrix(:,iObjNested)=sort(unique(gridOriginal(1:end,krigingIndexNested)),'descend');
            end
        end
    end
% -------------------------------------------------------------------------
    function  [gridToLow,gridToBig]=determinConfidenceBounds()
        % Define wanted quantile levels
        cumulativeProbability = @(x)(1/2)*(1+erf(x/sqrt(2)));
        enclosedProbalitiy = cumulativeProbability(obj.WidthConfidenceInterval) - cumulativeProbability(-obj.WidthConfidenceInterval);
        quantileValueLB = 0.5-enclosedProbalitiy/2;
        quantileValueUB = 0.5+enclosedProbalitiy/2;
        
        % Choose all points outside off and at the broder off desired
        % confidence region
        qToLow = pCover<=quantileValueLB;
        qToBig = pCover>=quantileValueUB;
        
        % Define confidence borders
        gridToLow = determineParetoSet_Mex(bsxfun(@times,gridOriginal(qToLow,:),obj.MinMax(krigingObjVector)));
        gridToLow = bsxfun(@times,gridToLow,obj.MinMax(krigingObjVector));
        gridToBig = determineParetoSet_Mex(bsxfun(@times,gridOriginal(qToBig,:),-obj.MinMax(krigingObjVector)));
        gridToBig = bsxfun(@times,gridToBig,-obj.MinMax(krigingObjVector));
        
    end
% -------------------------------------------------------------------------
    function  []=labelAxis(krigingObjVector)
        xlabel(obj.KrigingObjectNames{krigingObjVector(1)},'FontSize',20);
        ylabel(obj.KrigingObjectNames{krigingObjVector(2)},'FontSize',20);
        if nObj==3
            zlabel(obj.KrigingObjectNames{krigingObjVector(3)},'FontSize',20);
        end    
    end
    function [] = addLegend()
        warning('off','MATLAB:legend:IgnoringExtraEntries')
        hObjForLegend = [];
        cellStrLegend = {};
        for iObjNested = 1:6
            switch iObjNested
                case 1
                    if ~isempty(hExpected)
                        hObjForLegend(end+1) = hExpected;
                        cellStrLegend{end+1} = 'Expected Pareto Curve';
                    end
                case 2
                    if ~isempty(hConfidence)
                        hObjForLegend(end+1) = hConfidence;
                        cellStrLegend{end+1} = 'Confidence Tube';
                    end
                case 3
                    if ~isempty(objectOpt)
                        hObjForLegend(end+1) = objectOpt;
                        cellStrLegend{end+1} = 'Optimal Sample Points';
                    end
                case 4 
                    if ~isempty(objectNonOpt)
                        hObjForLegend(end+1) = objectNonOpt;
                        cellStrLegend{end+1} = 'Non-Optimal Sample Points';
                    end
                case 5
                    if ~isempty(objectOptExtra)
                        hObjForLegend(end+1) = objectOptExtra;
                        cellStrLegend{end+1} = 'Optimal Added Points';
                    end
                case 6
                    if ~isempty(objectNonOptExtra)
                        hObjForLegend(end+1) = objectNonOptExtra;
                        cellStrLegend{end+1} = 'Non-optimal Extra Points';
                    end
                otherwise
                    error('unknown graphical object')
            end
        end
        hLegend = legend(hObjForLegend,cellStrLegend);
        set(hLegend.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[.5;.5;.5;.5])); 
        warning('on','MATLAB:legend:IgnoringExtraEntries')
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
