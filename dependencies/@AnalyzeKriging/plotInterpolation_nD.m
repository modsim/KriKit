function [] = plotInterpolation_nD(obj,varargin)
% plotInterpolation_nD(KrigingObjectIndex)
%
% Uses the data which are calculated by calcInterpolation_nD and
% create a nD-Plot out of them
%
% Input: -
% - KrigingObjectIndex ... index of the kriging object which should
%                        be used in this function.
%
% Output: -
%
% You can set: -
%
% You can get: -
%
% For further details about settings see documentation of
% "plotScreeningAnalysis()"
% 
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

%% Initialization
KrigingObjectIndex = varargin{1};
KrigingObjectIndex = obj.checkKrigingIndizes(KrigingObjectIndex);
KrigingObjectIndexArray = KrigingObjectIndex;
KrigingObjectIndex = KrigingObjectIndexArray(1);

% Collecting Kriging predicion
Estimation = obj.KrigingPrediction_InterpolationnD{KrigingObjectIndex,1}(:,1);
minEstimation = min(Estimation);
maxEstimation = max(Estimation);

% Matrix containing all combination of input variables calues
% which were used for prediciton
InputMatrix = obj.KrigingPrediction_InterpolationnD{KrigingObjectIndex,2};
inputVarIndices = obj.KrigingPrediction_InterpolationnD{KrigingObjectIndex,3};

% Cell array containing the unique values of the input
% variables
InputMatrixProto = obj.KrigingPrediction_InterpolationnD{KrigingObjectIndex,4};
nRows = size(InputMatrixProto,2);

offsetIndices = 0;
if obj.getShowColorBar
    offsetIndices = 1;
end
%% Actual Plotting
figure();
hold on


% Plot Contour plots
for iPlots1=1:nRows
    for iPlots2=1:obj.getnPlots
        createContourPlot(iPlots1,iPlots2);
        set(gca,'PlotBoxAspectRatio',[1,1,1])
    end
end

%% Final Things

% Label plots 
for iPlots1=1:nRows
    labelPlotyAxis(iPlots1);
end
for iPlots1=1:nRows
    for iPlots2=1:obj.getnPlots
        labelPlotxAxis(iPlots1,iPlots2);
    end
end

% Title
if obj.getShowColorBar||obj.getnPlots>1
    subplot(nRows,obj.getnPlots+offsetIndices,ceil(obj.getnPlots/2))
end
title(obj.getKrigingObjectNames{KrigingObjectIndex},'FontSize',20)

% Colorbar
if obj.getNormColors==1
    caxis([minEstimation,maxEstimation])
end
if obj.getShowColorBar==1&&obj.getNormColors==1
    subplot(obj.getnPlots,obj.getnPlots+offsetIndices,obj.getnPlots+offsetIndices:obj.getnPlots+offsetIndices:obj.getnPlots*(obj.getnPlots+offsetIndices))
    colorbar
    axis off
end

% Adjust Color map
colormap(gcf,obj.ColormapToolbox);

%% ----------------------Nested Function -----------------------
% labelPlotxAxis
function [] = createContourPlot(iRow,iPlots2)

% Indices where the input values equal to the expected ones
if size(inputVarIndices,2)>2
    indexEstimation = InputMatrix(obj.Accuracy^2*obj.nPlots*(iRow-1)+1:obj.Accuracy^2*obj.nPlots*(iRow),...
                      inputVarIndices(iRow,3))==InputMatrixProto{inputVarIndices(iRow,3),iRow}(iPlots2);
else
    indexEstimation = true(size(InputMatrix,1),1);
end


% Choose Subplot window
indexSubPlot = (iRow-1)*(obj.nPlots+offsetIndices)+iPlots2;
subplot(nRows,obj.nPlots+offsetIndices,indexSubPlot);
hold on

% Actual Plotting
OutputPart = Estimation(obj.Accuracy^2*obj.nPlots*(iRow-1)+1:obj.Accuracy^2*obj.nPlots*(iRow));
OutputPart = OutputPart(logical(indexEstimation));
if inputVarIndices(1)<inputVarIndices(2)
    contourf(InputMatrixProto{inputVarIndices(iRow,1),iRow}',InputMatrixProto{inputVarIndices(iRow,2),iRow}',reshape(OutputPart,obj.Accuracy,obj.Accuracy)',obj.nContourLevels)
else
    contourf(InputMatrixProto{inputVarIndices(iRow,1),iRow}',InputMatrixProto{inputVarIndices(iRow,2),iRow}',reshape(OutputPart,obj.Accuracy,obj.Accuracy)',obj.nContourLevels)
end

% Adjusting color range
if obj.NormColors==1
    caxis([minEstimation,maxEstimation])
elseif obj.ShowColorBar==1
    colorbar
end

% Plot provided Data
if obj.ShowData
    for iKriging=KrigingObjectIndexArray'
        Data = obj.KrigingObjects{iKriging}.getInputData;
        plotSamples(Data,iRow,iPlots2,[255/255 102/255 0/255])
    end
end

end
% -------------------------------------------------------------
function []=plotSamples(Data,iRow,iPlots2,colorSamples)
% Plot only data which actual meet the desired values

% Copy Data
Data2 = Data;

% Choose the descrete varied input variable. Do not consider the
% continously varied variable (they comprises the entire measured range)
InputData = InputMatrixProto;
InputData = InputData{inputVarIndices(iRow,1:2),iRow}; % The remaining variable have only one value in each row

% Find Rows
r = (true(size(Data2,1),1));

% If length is smaller than 3 -> only two variables are available and
% therefore no one has to be set to a const value
if length(inputVarIndices)>2
    r=r&Data2(:,inputVarIndices(iRow,3))==InputMatrixProto{inputVarIndices(iRow,3),iRow}(iPlots2);
end
for iInputVar = setdiff(1:obj.KrigingObjects{KrigingObjectIndex}.getnInputVar,inputVarIndices(iRow,:));
    r=r&Data2(:,iInputVar)==InputMatrixProto{iInputVar,iRow};
end

% Actual Plotting
indicesinputVar = [obj.KrigingPrediction_InterpolationnD{KrigingObjectIndex,3}(iRow,1),...
                   obj.KrigingPrediction_InterpolationnD{KrigingObjectIndex,3}(iRow,2)];
plot(Data(r,indicesinputVar(1)),...
      Data(r,indicesinputVar(2)),...
      'ko','MarkerFaceColor',colorSamples);
axis([obj.LBInputVarInterpolation{KrigingObjectIndex}(indicesinputVar(1)),...
      obj.UBInputVarInterpolation{KrigingObjectIndex}(indicesinputVar(1)),...
      obj.LBInputVarInterpolation{KrigingObjectIndex}(indicesinputVar(2)),...
      obj.UBInputVarInterpolation{KrigingObjectIndex}(indicesinputVar(2))])
end
% -------------------------------------------------------------
function []=labelPlotxAxis(iRow,iPlots2)
    indexSubPlot = (iRow-1)*(obj.nPlots+offsetIndices)+iPlots2;
    subplot(nRows,obj.nPlots+offsetIndices,indexSubPlot)
    if isempty(obj.InputVarNames{KrigingObjectIndex})
        str1 = sprintf('InputVar %i',inputVarIndices(iRow,1));
        if size(inputVarIndices,2)>2
            str2=sprintf('\nInputVar %i = %0.2g',inputVarIndices(iRow,3),InputMatrixProto{inputVarIndices(iRow,3),iRow}(iPlots2));
        else 
            str2='';
        end
    else
        str1 = obj.InputVarNames{KrigingObjectIndex}(inputVarIndices(iRow,1));
        str1 = sprintf('\n\t\t%s',str1{1});
        if size(inputVarIndices,2)>2
            if size(inputVarIndices,2)>2
                str2=obj.InputVarNames{KrigingObjectIndex}(inputVarIndices(iRow,3));
            end
            str2=sprintf('\n%s = %0.2g',str2{1},InputMatrixProto{inputVarIndices(iRow,3),iRow}(iPlots2));
        else 
            str2 = '';
        end
    end
    xlabel(strcat(str1,str2),'FontSize',5*3/obj.nPlots+9);
end
% -------------------------------------------------------------
function []=labelPlotyAxis(iRow)

    indexSubPlot = (obj.nPlots+offsetIndices)*(iRow-1)+1;
    subplot(nRows,obj.nPlots+offsetIndices,indexSubPlot)
    if isempty(obj.InputVarNames{KrigingObjectIndex})
        str = sprintf('InputVar %i\n',inputVarIndices(iRow,2));
    else
        str = obj.InputVarNames{KrigingObjectIndex}(inputVarIndices(iRow,2));
        str = sprintf('%s',str{1});
    end
    ylabel(str,'FontSize',5*3/nRows+9);
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
