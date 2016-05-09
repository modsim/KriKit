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
% For triangular plotting the Ternplot package was used available at 
% http://www.mathworks.com/matlabcentral/fileexchange/2299-ternploton
% Acces Date (02-02-2016)
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
% Label plots on the left edge
for iPlots1=1:nRows
    for iPlots2=1:obj.getnPlots
        labelPlotxAxis(iPlots1,iPlots2);
    end
end

% Title
subplot(nRows,obj.getnPlots+1,ceil(obj.getnPlots/2))
title(obj.getKrigingObjectNames{KrigingObjectIndex},'FontSize',20)

% Colorbar
subplot(obj.getnPlots,obj.getnPlots+1,obj.getnPlots+1:obj.getnPlots+1:obj.getnPlots*(obj.getnPlots+1))
if obj.getNormColors==1
    caxis([minEstimation,maxEstimation])
end
if obj.getShowColorBar==1&&obj.getNormColors==1
    colorbar
end

% Adjust Color map
colormap(gcf,obj.ColormapToolbox);
axis off


%% ----------------------Nested Function -----------------------
% labelPlotxAxis
function [] = createContourPlot(iRow,iPlots2)
    % Initialization 
    startIndex = obj.getAccuracy^2*obj.getnPlots*(iRow-1) + (iPlots2-1)*obj.getAccuracy^2 + 1;
    endIndex = obj.getAccuracy^2*obj.getnPlots*(iRow-1) + iPlots2*obj.getAccuracy^2;
    validIndicesSub = startIndex:endIndex;
    
    % Indices where the input values are equal to the expected ones
    indexSubPlot = (iRow-1)*(obj.getnPlots+1)+iPlots2;
    
    % Plot
    subplot(nRows,obj.getnPlots+1,indexSubPlot);
    hold on
    [uniqueRow,indexRows]=unique(InputMatrix(validIndicesSub,:),'rows');
    ternpcolor(uniqueRow(:,inputVarIndices(iRow,1)),...
               uniqueRow(:,inputVarIndices(iRow,2)),...
               uniqueRow(:,inputVarIndices(iRow,3)),...
                Estimation(validIndicesSub(indexRows)))
    if isempty(obj.getInputVarNames(KrigingObjectIndex(1)))
        xString = sprintf('Var %i=%g\n',inputVarIndices(iRow,1),max(uniqueRow(:,inputVarIndices(iRow,1))) );
        yString = sprintf('\n\n\n\t\t                       Var %i=%g',inputVarIndices(iRow,2),max(uniqueRow(:,inputVarIndices(iRow,2))) );
        zString = sprintf('\n\n\nVar %i=%g\t\t                       ',inputVarIndices(iRow,3),max(uniqueRow(:,inputVarIndices(iRow,3))) );

    else
        stringInputVar = obj.getInputVarNames(KrigingObjectIndex(1));
        xString = sprintf('\n\n\n%s=%g\n',stringInputVar{inputVarIndices(iRow,1)},max(uniqueRow(:,inputVarIndices(iRow,1))) );
        yString = sprintf('%s=%g',stringInputVar{inputVarIndices(iRow,2)},max(uniqueRow(:,inputVarIndices(iRow,2))) );
        zString = sprintf('\n\n\n\t\t                       %s=%g\t\t                       ',stringInputVar{inputVarIndices(iRow,3)},max(uniqueRow(:,inputVarIndices(iRow,3))) );
    end
    ternlabel(zString,yString,xString);
    
    if obj.getNormColors==1
        caxis([minEstimation,maxEstimation])
    elseif obj.getShowColorBar==1
        colorbar
    end

    % Plot provided Data
    if obj.getShowData==1
        for iKriging=KrigingObjectIndexArray'
            Data = obj.KrigingObjects{iKriging}.getInputData;
            plotSamples(Data,iRow,iPlots2)
        end
    end

end
% -------------------------------------------------------------
function []=plotSamples(Data,iRow,iPlots2)
%                 Data = obj.KrigingObjects{KrigingObjectIndex}.getInputData;
    % Plot only data which actual met the desired values
        % Copy Data
    Data2 = Data;
    % Choose the descrete varied input variable

    % Do not consider the continously varied
    % variable (they comprises the entire measured range)
    InputData = InputMatrixProto;
    InputData = InputData{inputVarIndices(iRow,1:2),iRow}; % The remaining variable have only one value in each row
    
        % Find Rows
    r = (true(size(Data2,1),1));

    % If length is smaller than 3 -> only two variables are avaiable and
    % therefore no one has to be set on a const value
    if length(inputVarIndices)>3
        r=r&Data2(:,inputVarIndices(iRow,4))==InputMatrixProto{inputVarIndices(iRow,4),iRow}(iPlots2);
    end

    for iInputVar = setdiff(1:obj.KrigingObjects{KrigingObjectIndex}.getnInputVar,inputVarIndices(iRow,:));
        r=r&Data2(:,iInputVar)==InputMatrixProto{iInputVar,iRow};
    end

   ternplot(Data(r,obj.KrigingPrediction_InterpolationnD{KrigingObjectIndex,3}(iRow,1)),...
            Data(r,obj.KrigingPrediction_InterpolationnD{KrigingObjectIndex,3}(iRow,2)),...
            Data(r,obj.KrigingPrediction_InterpolationnD{KrigingObjectIndex,3}(iRow,3)),...
            'ko','MarkerFaceColor',[255/255 102/255 0/255]);
end
% -------------------------------------------------------------
function []=labelPlotxAxis(iRow,iPlots2)
%     subplot(obj.getnPlots,obj.getnPlots+1,(obj.getnPlots-1)*(obj.getnPlots+1)+iPlots1)
    indexSubPlot = (iRow-1)*(obj.getnPlots+1)+iPlots2;
    subplot(nRows,obj.getnPlots+1,indexSubPlot)
    
    % Use default names of input variables when user did not provide any
    if isempty(obj.getInputVarNames(KrigingObjectIndex))
        str1 = sprintf('\n');
        
        % Create Label
        if size(inputVarIndices,2)>2
            str2=sprintf('\nInputVar %i = %0.2g',inputVarIndices(iRow,4),InputMatrixProto{inputVarIndices(iRow,4),iRow}(iPlots2));
        else 
            str2='';
        end
    else
        inputNames = obj.getInputVarNames(KrigingObjectIndex);
        str1 = sprintf('\n\n\t\t%s');
        if size(inputVarIndices,2)>2
            if size(inputVarIndices,2)>2
                str2=inputNames(inputVarIndices(iRow,4));
            end
            str2=sprintf('\n%s = %0.2g',str2{1},InputMatrixProto{inputVarIndices(iRow,4),iRow}(iPlots2));
        else 
            str2 = '';
        end
    end
    xlabel(strcat(str1,str2),'FontSize',5*3/obj.getnPlots+9);
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
