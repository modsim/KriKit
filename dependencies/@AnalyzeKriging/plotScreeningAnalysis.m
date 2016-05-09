function [] = plotScreeningAnalysis(obj,varargin)
% [] = plotScreeningAnalysis(obj,KrigingObjectIndex,Objective)
% KrigingObjectIndex ... index of the Kriging object which should
%                        be used here
% This function uses the data generated using "calcScreeningAnalysis" and
% generates a subplot containing the contour plot of all possible
% combinations of the input variables
%
% % Input:
% - KrigingObjectIndex ... indices of kriging objects of interest [1XnObj]
%
% - Objective ... decide if the Kriging estimation ('KrigingInterpolation'), the expected
%                 improvement ('ExpectedImprovement') or optimal regions
%                 ('Optimum') shall be plotted
%
% 
% Output: -
% 
% You can set:
% nContourLevels ... Defines how many levels are drawn in the subplots. By
%                    default 10.
% ShowColorBar   ... If true colorbar is shown. By default true
%
% NormColors  ... If true, color of countours in all subplots are
%                 associated with the same values. The range of the
%                 contours is defined by [Min(Prediction),Max(Prediction)].
%                 If false, each subplot has its own color bar,
% You can get: -
%
%
% For further details about settings w.r.t. expected improvement or
% identification of optimal regions see documentation of
% "calcExpectedImprovementFromPredictions", or "plotOptimum23D()",
% respectively
% 
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.


%% Initialization
KrigingObjectIndex = varargin{1};
if length(KrigingObjectIndex)~=1
    error('KrigingObjectIndex must be a scalar')
end
Objective = varargin{2};
[nCombinations,nPossibleCombi,minEstimation,maxEstimation,plottingObjective] = doInitialiaztion(KrigingObjectIndex,Objective);

%% Actual Plotting
figure()
if nPossibleCombi>=3
    subplot(nPossibleCombi,nPossibleCombi,floor((nPossibleCombi+1)/2):ceil((nPossibleCombi+1)/2))
else
    subplot(nPossibleCombi,nPossibleCombi,ceil((nPossibleCombi+1)/2):ceil((nPossibleCombi+1)/2))
end
axis off
title(obj.KrigingObjectNames{KrigingObjectIndex},'FontSize',20)
for iCombination = 1 : nCombinations
    createContourPlot(iCombination)
end


%% Final Adjustements
labelPlotAxis

nInputVar=obj.KrigingObjects{KrigingObjectIndex}.getnInputVar;
if obj.NormColors==1&&obj.ShowColorBar==1&&nInputVar>2
    subplot(nPossibleCombi,nPossibleCombi,nInputVar-1:nPossibleCombi:(nInputVar-2-1)*nPossibleCombi+obj.KrigingObjects{KrigingObjectIndex}.getnInputVar-1);
    axis off
    caxis([minEstimation,maxEstimation])
    colorbar('location','East')
elseif obj.NormColors==1&&obj.ShowColorBar==1&&nInputVar==2
    subplot(nPossibleCombi,nPossibleCombi,1);
    axis off
    caxis([minEstimation,maxEstimation])
    colorbar('location','East')
end

%% Nested Functions
function [nCombinations,nPossibleCombi,minEstimation,maxEstimation,plottingObjective] = doInitialiaztion(KrigingObjectIndex,Objective)
    nCombinations = length(obj.KrigingPrediction_Screening{KrigingObjectIndex,1});
    nPossibleCombi = (obj.KrigingObjects{KrigingObjectIndex}.getnInputVar-1);
    minEstimation = inf;
    maxEstimation = -inf;
    for iCombinationNested = 1 : nCombinations
        switch Objective
            case 'KrigingInterpolation'
                minEstimation = min(minEstimation,min(obj.KrigingPrediction_Screening{KrigingObjectIndex,1}{iCombinationNested}(:,1)));
                maxEstimation = max(maxEstimation,max(obj.KrigingPrediction_Screening{KrigingObjectIndex,1}{iCombinationNested}(:,1)));
                plottingObjective = [];
            case 'ExpectedImprovement'
                plottingObjective = obj.calcExpectedImprovementFromPredictions(KrigingObjectIndex,obj.KrigingPrediction_Screening{KrigingObjectIndex,1}{iCombinationNested});
                minEstimation = min(minEstimation,min(plottingObjective));
                maxEstimation = max(maxEstimation,max(plottingObjective));
            case 'Optimum'
                minEstimation = 0;
                maxEstimation = 1;
                plottingObjective = [];
            otherwise
                error('Unknown plotting objective')
        end
    end
end
% -------------------------------------------------------------
function [] = createContourPlot(iCombination)
    switch Objective
        case 'KrigingInterpolation'
            plottingObjective = obj.KrigingPrediction_Screening{KrigingObjectIndex,1}{iCombination};
        case 'ExpectedImprovement'
            prediction = obj.KrigingPrediction_Screening{KrigingObjectIndex,1}{iCombination};
            prediction(:,1) = -obj.MinMax(KrigingObjectIndex)*prediction(:,1);
            plottingObjective = obj.calcExpectedImprovementMainPart(KrigingObjectIndex,prediction);

        case 'Optimum'
            if length(varargin)>2
                testValue = varargin{3};
            else
                testValue = [];
            end
            [inputDataOpt,indexValid] = doTestForOptimality(obj,KrigingObjectIndex,4,testValue);
             indexInvalid = ~indexValid;  

        otherwise
            error('Unknown plotting objective')
    end
    
    % Indices of input variables plotted on the x- and y-axis in the
    % current subplot
    indexInputVar = obj.KrigingPrediction_Screening{KrigingObjectIndex,3}{iCombination};
    indexSubPlot = (indexInputVar(2)-2)*nPossibleCombi + indexInputVar(1);

    % Actual plotting
    subplot(nPossibleCombi,nPossibleCombi,indexSubPlot);
    hold on
    switch Objective
        case {'KrigingInterpolation','ExpectedImprovement'}
            inputData = obj.KrigingPrediction_Screening{KrigingObjectIndex,2}{iCombination};
            contourf(unique(inputData(:,indexInputVar(1))),unique(inputData(:,indexInputVar(2))),...
                reshape(plottingObjective(:,1),obj.Accuracy,obj.Accuracy)',obj.nContourLevels)
            
            % Set Colormap
            colormap(gcf,obj.ColormapToolbox);
        case 'Optimum'
            % Fill in the t-test results
            hold on
            plot(inputDataOpt(indexValid,1),inputDataOpt(indexValid,2),'r*');
            plot(inputDataOpt(indexInvalid,1),inputDataOpt(indexInvalid,2),'b*');
            axis tight
        otherwise
        error('Unknown plotting objective')
    end
    
    if obj.ShowData
        plotDataAndMarkOutliers(4,KrigingObjectIndex,iCombination)
    end
            
    if obj.NormColors==1
        caxis([minEstimation,maxEstimation])
    elseif obj.ShowColorBar==1
        colorbar
    end
end
% -------------------------------------------------------------
function []=labelPlotAxis()
    % xlabels
    for iInputVar = 1 : obj.KrigingObjects{KrigingObjectIndex}.getnInputVar-1
        % Last row
        subplot(nPossibleCombi,nPossibleCombi,(obj.KrigingObjects{KrigingObjectIndex}.getnInputVar-2)*nPossibleCombi+iInputVar)
        if isempty(obj.InputVarNames{KrigingObjectIndex})
%             xlabel(sprintf('InputVar %i',inputVarIndices(2)))
            xlabel(sprintf('InputVar %i',iInputVar))
        else
            xlabel(obj.InputVarNames{KrigingObjectIndex}(iInputVar))
        end
        set(gca,'xaxislocation','bottom');
    end
    
    % ylabels
    for iInputVar = 2 : obj.KrigingObjects{KrigingObjectIndex}.getnInputVar
        % First Column
        subplot(nPossibleCombi,nPossibleCombi,(iInputVar-2)*nPossibleCombi + 1)
        if isempty(obj.InputVarNames{KrigingObjectIndex})
%             ylabel(sprintf('InputVar %i',inputVarIndices(2)))
            ylabel(sprintf('InputVar %i',iInputVar))
        else
            ylabel(obj.InputVarNames{KrigingObjectIndex}(iInputVar))
        end
        set(gca,'yaxislocation','left');
    end
    

end
% ------------------------------------------------------------------------
function [] = plotDataAndMarkOutliers(dimensionInterpolation,KrigingObjectIndexArray,iCombination)
    for iKrigingObject=KrigingObjectIndexArray'
        [Data,dataShown,OutlierShown]=obj.determineRelevantDataPointsForPlot(iKrigingObject,dimensionInterpolation,iCombination);
        
        indicesInputVar = obj.KrigingPrediction_Screening{iKrigingObject,3}{iCombination};
        plot(Data(dataShown,indicesInputVar(1)),...
        Data(dataShown,indicesInputVar(2)),...
        'ko','MarkerFaceColor',[255/255 102/255 0/255]);

        plot(Data(OutlierShown,indicesInputVar(1)),...
             Data(OutlierShown,indicesInputVar(2)),...
        'ko','MarkerFaceColor','k');

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
