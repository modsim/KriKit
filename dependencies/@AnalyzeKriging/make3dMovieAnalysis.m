function [] = make3dMovieAnalysis(obj,varargin)
% [] = make3dMovieAnalysis(obj,KrigingObjectIndex,NonContiousInputVar,VaryingInputVar,FileName)
%
% Input:
% - KrigingObjectIndex ... Index of the KrigingObject of interest
%
% - NonContiousInputVar ...  (Array) containing All indices of input
%                            variables which are not plotted on the x- or
%                            y-axis of the plot 
%                            (1X(nInputVar-number of variables continously varied on the axis))
% - VaryingInputVar     ...  (Array) containing indices of input variables
%                            which are changing during the video.
%                            VaryingInputVar is a subset of
%                            NonContiousInputVar 
% - FileName            ... (String) Name of the video-file e.g. 'test.avi'   
%
% Output: -
%
% You can set(optional):
%
% - ReferencePoint      ... (Array) Contains the fixed  values of the input
%                           variables defined in "NonContiousInputVar". If
%                           not set then the parameter set with the highest
%                           number of samples is chosen (for more details
%                           see documentation of
%                           "getUniqueInputVarCombinations()")
% - FrameRate           ... (int) Number of images in one second
% - nStepsBetweenBounds ... (int) Intermediate steps between the bound given 
%                           in "LB(UB)InterpolationRange". By default 1(for
%                           more details about "LB(UB)InterpolationRange"
%                           see documentation of
%                           "calcMutualInterpolation_23D()")
% - PlottingRange       ... (Array) Ranges of the y.axis(2D-plot) or
%                           z-axis(3D plot). By default
%                           [min(outputData),max(outputData)]
%
% You can get:-
%
% For further otion regarding interpolation and plotting setting see
% documentation of "calcMutualInterpolation_23D()" and ""
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

% Initialization
KrigingObjectIndex = varargin{1};
NonContiousInputVar = varargin{2};
VaryingInputVar = varargin{3};
FileName = varargin{4};

backUpSuppressFigure=obj.SuppressFigure;
[continousInputVar,ReferencePoint,range,plottingRange] = doInitialization(NonContiousInputVar,VaryingInputVar);
dimensionPlot = length(continousInputVar)+1;
if dimensionPlot~=2&&dimensionPlot~=3
    error('nInputvar - length(NonContiousInputVar) has to be either 1(2D plot) or 2(3D plot), but it is %i',length(continousInputVar))
end

% Create InputVar-Combinations for the several images used for the movie
[InputValues]=createInputVariableCombinations(VaryingInputVar);

% Create the movie
[writerObj]=createMovie();
% Final stuff
close(writerObj);
delete('tmp.png')
obj.SuppressFigure=backUpSuppressFigure;
%% Nested Functions
    function [continousInputVar,ReferencePoint,range,plottingRange] = doInitialization(NonContiousInputVar,VaryingInputVar)
        if isempty(NonContiousInputVar)
            error('"NonContiousInputVar" is not defined. Use its setNonContiousInputVar')
        end

        continousInputVar=setdiff(1: obj.KrigingObjects{1}.getnInputVar,NonContiousInputVar);
%         continousInputVar = obj.getVaryingInputVar;

        if isempty(continousInputVar)
            error('"NonContiousInputVar" contains to many indices. max(length(NonContiousInputVar))=nInputvar-1')
        end

        % The point of the var which are changing during the
        % video are irrelevant since the given range are
        % concerned
        if isempty(obj.ReferencePoint)
            obj.calcAndPlotInterpolation_3D_BestChoice(KrigingObjectIndex,continousInputVar);
            ReferencePoint = [obj.getChosenCombinationsForPlot];
            close all
        else
            if length(obj.ReferencePoint)~=length(NonContiousInputVar)&&length(NonContiousInputVar)>1
                error('Length of "ReferencePoint" and "NonContiousInputVar" are not the same')
            end
            ReferencePoint = obj.ReferencePoint;   % contains only reference values of the noncontinousvar. 
        end
        % Input variable which is changing during the video 
        if isempty(VaryingInputVar); 
            error('"VaryingInputVar" is not defined. Use its setVaryingInputVar')
        end
        
        defineBoundOfInputVar(obj,KrigingObjectIndex);
        range = [obj.LBInputVarInterpolation{KrigingObjectIndex}',...
                 obj.UBInputVarInterpolation{KrigingObjectIndex}'];


        if isempty(FileName)||~ischar(FileName)
            error('"FileName" is not defined. Use its setFileName')
        end

        if isempty(obj.PlottingRange)
            plottingRange = ([min(obj.KrigingObjects{KrigingObjectIndex}.getOutputData),...
                       max(obj.KrigingObjects{KrigingObjectIndex}.getOutputData)]);
        else
            plottingRange = obj.PlottingRange;
        end
    end
% -------------------------------------------------------------------------
    function [InputValues]=createInputVariableCombinations(VaryingInputVar)
        % Initialization
        InputValuesCell = cell(length(NonContiousInputVar),1);
        InputValues = ones(obj.nStepsBetweenBounds,length(NonContiousInputVar));
        
        % Decide in which range the input variables should be varied
        rangeProto = [min(obj.KrigingObjects{1}.getInputData);max(obj.KrigingObjects{1}.getInputData)]';
        rangeProto(VaryingInputVar,:)=range(VaryingInputVar,:);
        rangeProto(continousInputVar,:)=[];

        % Create grid for each input variable
        for iInputVar = 1: length(NonContiousInputVar)
            if sum(VaryingInputVar==NonContiousInputVar(iInputVar))>0
                InputValuesCell{iInputVar} = linspace(rangeProto(iInputVar,1),rangeProto(iInputVar,2),obj.nStepsBetweenBounds);
            else
                InputValuesCell{iInputVar} = ReferencePoint(iInputVar);
            end
        end

        % Combine grids
        for iInputVar = 1: length(NonContiousInputVar)
            InputValues(:,iInputVar) = InputValues(:,iInputVar).*InputValuesCell{iInputVar}(:);
        end
    end
% -------------------------------------------------------------------------
    function [writerObj]=createMovie()
        for iStep = 1:obj.nStepsBetweenBounds
            % Provide information about progress
            fprintf('Step %i of %i\n',iStep,obj.nStepsBetweenBounds)
            
            % Calculate the interpolation
            obj.calcMutualInterpolation_23D(KrigingObjectIndex,continousInputVar,NonContiousInputVar,InputValues(iStep,:),dimensionPlot);
            obj.plotInterpolation_23D(KrigingObjectIndex,dimensionPlot)
            grid on
            
            h = adjustPlot(iStep);
            
            % Capture figure 
            if iStep == 1
                if dimensionPlot>2
                    obj.SuppressFigure=0;
                end
                [cameraPosition,writerObj]=doInitialCapture(h);
            else
                obj.SuppressFigure=1;
                if dimensionPlot>2
                    campos(cameraPosition);
                end
                saveas(h,'tmp.png','png');
                frame = im2frame(imread('tmp.png'));
                writeVideo(writerObj,frame);
            end
            close(gcf)
        end
    end
% -------------------------------------------------------------------------
    function [h]=adjustPlot(iStep)
        if(length(VaryingInputVar)==1)
            title(...
                sprintf('%s: %9.2e',obj.InputVarNames{KrigingObjectIndex(1)}{VaryingInputVar},InputValues(iStep,VaryingInputVar==NonContiousInputVar)),...
                'FontSize',obj.FontSize)
        end
        if dimensionPlot==3
            zlim(plottingRange)
        elseif dimensionPlot==2
            ylim(plottingRange)
        else 
            error('Unexpected plot dimension%g',dimensionPlot)
        end
        h = gcf;
    end
% -------------------------------------------------------------------------
    function [cameraPosition,writerObj] = doInitialCapture(h)
        % Set properties
        set(gca,'nextplot','replacechildren');
        set(gcf,'Renderer','zbuffer');
        
        % Let user decide about the perspective of the 3D-plot
        if dimensionPlot>2
            disp('Please Rotate Figure and continue with an arbitrary key')
            pause
        end
        cameraPosition=campos();
        
        % Make picture from curent plot
        saveas(h,'tmp.png','png');
        
        % open/create new file and save the picture
        writerObj = VideoWriter(FileName);
        set(writerObj,'FrameRate',obj.FrameRate)
        open(writerObj);
        frame = im2frame(imread('tmp.png'));
        writeVideo(writerObj,frame);
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
