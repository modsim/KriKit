function varargout = KrigingGUI(varargin)
% KRIGINGGUI MATLAB code for KrigingGUI.fig
%      KRIGINGGUI, by itself, creates a new KRIGINGGUI or raises the existing
%      singleton*.
%
%      H = KRIGINGGUI returns the handle to a new KRIGINGGUI or the handle to
%      the existing singleton*.
%
%      KRIGINGGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in KRIGINGGUI.M with the given input arguments.
%
%      KRIGINGGUI('Property','Value',...) creates a new KRIGINGGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before KrigingGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to KrigingGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

% Edit the above text to modify the response to help KrigingGUI

% Last Modified by GUIDE v2.5 17-Jun-2015 14:53:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @KrigingGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @KrigingGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before KrigingGUI is made visible.
function KrigingGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to KrigingGUI (see VARARGIN)
if exist('KrigingAnalysisObj')~=0
    clear KrigingAnalysisObj
end

% Overall Kriging ANalysis Object
handles.KrigingAnalysisObj = AnalyzeKriging;

% Choose default command line output for KrigingGUI
handles.output = hObject;



f = gcf;
a = axes;
set(a, 'Visible', 'off');
% # Stretch the axes over the whole figure.
set(a, 'Position', [0, 0, 1, 1]);
% # Switch off autoscaling.
set(a, 'Xlim', [0, 1], 'YLim', [0, 1]);

% # Create all the controls.
% uicontrol('Parent', f, 'Style', 'edit', 'String', 'Input...');

% # Draw!
line([0.05, 0.95], [0.67, 0.67], 'Parent', a,'Color',[0.5,0.5,0.5],'LineWidth',1.5)
line([0.05, 0.95], [0.365, 0.365], 'Parent', a,'Color',[0.5,0.5,0.5],'LineWidth',1.5)



% Update handles structure
guidata(hObject, handles);

% UIWAIT makes KrigingGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = KrigingGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in AddObjective.
function AddObjective_Callback(hObject, eventdata, handles)
% hObject    handle to AddObjective (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
nameStr = inputdlg('Enter the name of your objectve','Add new Objective');

if ~isempty(nameStr)
    handles.KrigingAnalysisObj.addKrigingObject(2,nameStr);

    % Update List of Objectives
    nameStrProto = handles.KrigingAnalysisObj.getKrigingObjectNames;
    nameStr = cell(length(nameStrProto),1);
    for iObj=1:length(nameStrProto)
        nameStr{iObj} = nameStrProto{iObj}{1};
    end
    set(handles.KrigingObjectivePopUp,'String',nameStr);

    set(handles.KrigingObjectivePopUp,'Value',handles.KrigingAnalysisObj.getnKrigingObjects)
    set(handles.ListInputVar,'String',handles.KrigingAnalysisObj.getInputVarNames(get(handles.KrigingObjectivePopUp,'Value')))
    
    % Set Prefered conditions
    handles.KrigingAnalysisObj.KrigingObjects{length(nameStrProto)}.setUseInverse(true);
    handles.KrigingAnalysisObj.KrigingObjects{length(nameStrProto)}.setShowDetails(true);
    % get(handles.KrigingObjectivePopUp,'Value') = handles.KrigingAnalysisObj.getnKrigingObjects;
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in LoadDataButton.
function LoadDataButton_Callback(hObject, eventdata, handles)
% hObject    handle to LoadDataButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.KrigingAnalysisObj.KrigingObjects)
    errordlg('No Kriging Object chosen')
    return
end

inputData = [];
outputData = [];
[filenameIn,pathIn] = uigetfile('*.txt','Select Input Data');

if filenameIn~=0
    inputData = importdata(strcat(pathIn,filenameIn));
end

if ~isempty(inputData)
    [filenameOut,pathOut] = uigetfile('*.txt','Select Output Data');
    if filenameOut~=0
        outputData = importdata(strcat(pathOut,filenameOut));
    end
end

if ~isempty(inputData)&&~isempty(outputData)
    if isempty(handles.KrigingAnalysisObj.KrigingObjects{get(handles.KrigingObjectivePopUp,'Value')}.getInputData)
        handles.KrigingAnalysisObj.KrigingObjects{get(handles.KrigingObjectivePopUp,'Value')}.setNormInput(1)
    end
    if isempty(handles.KrigingAnalysisObj.KrigingObjects{get(handles.KrigingObjectivePopUp,'Value')}.getOutputData)
        handles.KrigingAnalysisObj.KrigingObjects{get(handles.KrigingObjectivePopUp,'Value')}.setNormOutput(1)
    end
    handles.KrigingAnalysisObj.KrigingObjects{get(handles.KrigingObjectivePopUp,'Value')}.setInputData(inputData)
    handles.KrigingAnalysisObj.KrigingObjects{get(handles.KrigingObjectivePopUp,'Value')}.setOutputData(outputData)
    handles.KrigingAnalysisObj.KrigingObjects{get(handles.KrigingObjectivePopUp,'Value')}.setCovariogramModelChoice(3);
    handles.KrigingAnalysisObj.KrigingObjects{get(handles.KrigingObjectivePopUp,'Value')}.setBasisFct('polynomial',0)
    if isempty(handles.KrigingAnalysisObj.getInputVarNames(get(handles.KrigingObjectivePopUp,'Value')))
        indices = 1:handles.KrigingAnalysisObj.KrigingObjects{get(handles.KrigingObjectivePopUp,'Value')}.getnInputVar;
        indicesCell = num2cell(indices);
        for iCell = 1:length(indicesCell)
            indicesCell{iCell} = num2str(indicesCell{iCell});
        end
        
%         indices = 1:handles.KrigingAnalysisObj.KrigingObjects{get(handles.KrigingObjectivePopUp,'Value')}.getnInputVar;
        handles.KrigingAnalysisObj.setInputVarNames(get(handles.KrigingObjectivePopUp,'Value'),indicesCell);
%         handles.KrigingAnalysisObj.setInputVarNames(get(handles.KrigingObjectivePopUp,'Value'),num2str(indices'));
    end
    set(handles.ListInputVar,'String',handles.KrigingAnalysisObj.getInputVarNames(get(handles.KrigingObjectivePopUp,'Value')))
%     handles.
end

% Update handles structure
guidata(hObject, handles);




% --- Executes on button press in DeleteObjective.
function DeleteObjective_Callback(hObject, eventdata, handles)
% hObject    handle to DeleteObjective (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str = inputdlg('Enter the index of your objectve you want to remove','Delete Objective');
if ~isempty(str)
    indexObj = str2double(str{1});

    if (~isempty(str)&&indexObj<=handles.KrigingAnalysisObj.getnKrigingObjects)
        handles.KrigingAnalysisObj.removeKrigingObject(indexObj);

        % Update List of Objectives
        nameStrProto = handles.KrigingAnalysisObj.getKrigingObjectNames;
        nameStr = cell(length(nameStrProto),1);
        for iObj=1:length(nameStrProto)
            nameStr{iObj} = nameStrProto{iObj}{1};
        end
        
        if isempty(nameStr)
            set(handles.KrigingObjectivePopUp,'String','Add Objective');
        else
            set(handles.KrigingObjectivePopUp,'String',nameStr);
        end

        set(handles.KrigingObjectivePopUp,'Value',1)
        if ~isempty(nameStr)
            set(handles.ListInputVar,'String',handles.KrigingAnalysisObj.getInputVarNames(get(handles.KrigingObjectivePopUp,'Value')))
        else
            set(handles.ListInputVar,'String',{})
        end
    end
end
% Update handles structure
guidata(hObject, handles);



% --- Executes on selection change in KrigingObjectivePopUp.
function KrigingObjectivePopUp_Callback(hObject, eventdata, handles)
% hObject    handle to KrigingObjectivePopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns KrigingObjectivePopUp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from KrigingObjectivePopUp

% Update handles structure
set(handles.ListInputVar,'String',handles.KrigingAnalysisObj.getInputVarNames(get(handles.KrigingObjectivePopUp,'Value')))
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function KrigingObjectivePopUp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to KrigingObjectivePopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in ListInputVar.
function ListInputVar_Callback(hObject, eventdata, handles)
% hObject    handle to ListInputVar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ListInputVar contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ListInputVar
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function ListInputVar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ListInputVar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% % --- Executes on button press in pushbutton4.
% function pushbutton4_Callback(hObject, eventdata, handles)
% % hObject    handle to pushbutton4 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in editInputVarNames.
function editInputVarNames_Callback(hObject, eventdata, handles)
% hObject    handle to editInputVarNames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nameStr = inputdlg('Type in New Name of The Input Variable','Change Name of Input Variable');
if ~isempty(nameStr)
    inputVarNames = handles.KrigingAnalysisObj.getInputVarNames(get(handles.KrigingObjectivePopUp,'Value'));
    inputVarNames{get(handles.ListInputVar,'Value')} = nameStr{1};
    handles.KrigingAnalysisObj.setInputVarNames(get(handles.KrigingObjectivePopUp,'Value'),inputVarNames);
    set(handles.ListInputVar,'String',handles.KrigingAnalysisObj.getInputVarNames(get(handles.KrigingObjectivePopUp,'Value')))
end

guidata(hObject, handles);


% --- Executes on button press in EstiCoVariogramPushButton.
function EstiCoVariogramPushButton_Callback(hObject, eventdata, handles)
% hObject    handle to EstiCoVariogramPushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


[handles.KrigingAnalysisObj] = dialogCovarEsti(handles);


% --------------------------------------------------------------------
function uipushtoolSaveFile_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtoolSaveFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uisave('handles')

% --------------------------------------------------------------------
function uipushtoolOpenFile_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtoolOpenFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% uiopen('*.mat')
[krigingFile,krigingPath] = uigetfile('*.mat');
if krigingFile~=0
    handlesOld=load(strcat(krigingPath,krigingFile));
    handles.KrigingAnalysisObj = handlesOld.handles.KrigingAnalysisObj;
    
    nameStrProto = handles.KrigingAnalysisObj.getKrigingObjectNames;
    nameStr = cell(length(nameStrProto),1);
    for iObj=1:length(nameStrProto)
        nameStr{iObj} = nameStrProto{iObj}{1};
    end
    set(handles.KrigingObjectivePopUp,'Value',1);
    set(handles.KrigingObjectivePopUp,'String',nameStr);
    
    set(handles.ListInputVar,'String',handles.KrigingAnalysisObj.getInputVarNames(get(handles.KrigingObjectivePopUp,'Value')))
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pushbutton3DInterpolation.
function pushbutton3DInterpolation_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3DInterpolation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dialogInterpolation(handles,3);


% --- Executes on button press in pushbutton2DInterpolation.
function pushbutton2DInterpolation_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2DInterpolation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dialogInterpolation(handles,2);


% --- Executes on button press in pushbuttonNDInterpolation.
function pushbuttonNDInterpolation_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonNDInterpolation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dialogInterpolation(handles,4);


% 
function [backUPShowDetails]=doInitialStuffForDataAnalysis(handles)
choice = questdlg('Would you like to save the results in a text file?','Save Results','Yes','No','No');

if strcmp(choice,'Yes');
    [file,path] = uiputfile('DataAnalysis.txt','Save file name');
end


% Save the state if infroamtion should be shown in the command window and
% temporarily enable this option
backUPShowDetails = handles.KrigingAnalysisObj.KrigingObjects{get(handles.KrigingObjectivePopUp,'Value')}.getShowDetails;
handles.KrigingAnalysisObj.KrigingObjects{get(handles.KrigingObjectivePopUp,'Value')}.setShowDetails(true);

% Overrite the file
if strcmp(choice,'Yes');
    fileID = fopen(strcat(path,file),'w');
    fclose(fileID);
    diary(strcat(path,file))
end
clc

% stringPath = strcat(path,file);

% --- Executes on button press in pushbuttonCCD.
function pushbuttonCCD_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCCD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Interpolation Type
% strNames = varargin{2};
if handles.KrigingAnalysisObj.getnKrigingObjects==0
    msgbox('No Kriging object was selected')
    error('No Kriging object was selected')
end

[backUPShowDetails]=doInitialStuffForDataAnalysis(handles);

% Check if Krikit or Matlab statistic toolbox shall be used
backUpGPR = handles.KrigingAnalysisObj.KrigingObjects{get(handles.KrigingObjectivePopUp,'Value')}.getUseMatlabRegressionGP;
handles.KrigingAnalysisObj.KrigingObjects{get(handles.KrigingObjectivePopUp,'Value')}.setUseMatlabRegressionGP(false);

handles.KrigingAnalysisObj.doCompositeDesignAnalysis(get(handles.KrigingObjectivePopUp,'Value'))

% Go back to the original state
handles.KrigingAnalysisObj.KrigingObjects{get(handles.KrigingObjectivePopUp,'Value')}.setUseMatlabRegressionGP(backUpGPR);


% Close it
diary('off')

% Reset current State
handles.KrigingAnalysisObj.KrigingObjects{get(handles.KrigingObjectivePopUp,'Value')}.setShowDetails(backUPShowDetails);


% --- Executes on button press in pushbuttonFFD.
function pushbuttonFFD_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonFFD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[backUPShowDetails]=doInitialStuffForDataAnalysis(handles);

% Check if Krikit or Matlab statistic toolbox shall be used
backUpGPR = handles.KrigingAnalysisObj.KrigingObjects{get(handles.KrigingObjectivePopUp,'Value')}.getUseMatlabRegressionGP;
handles.KrigingAnalysisObj.KrigingObjects{get(handles.KrigingObjectivePopUp,'Value')}.setUseMatlabRegressionGP(false);

handles.KrigingAnalysisObj.doFullFactorialAnalysis(get(handles.KrigingObjectivePopUp,'Value'))

% Go back to the original state
handles.KrigingAnalysisObj.KrigingObjects{get(handles.KrigingObjectivePopUp,'Value')}.setUseMatlabRegressionGP(backUpGPR);

% Close it
diary('off')

% Reset current State
handles.KrigingAnalysisObj.KrigingObjects{get(handles.KrigingObjectivePopUp,'Value')}.setShowDetails(backUPShowDetails);


% --- Executes on button press in pushbuttonScreening.
function pushbuttonScreening_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonScreening (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
outputRef = dialogReferencePoint(handles);
switch outputRef 
    case 'OK Button Was Pushed'
        currentObj = get(handles.KrigingObjectivePopUp,'Value');
        handles.KrigingAnalysisObj.calcScreeningAnalysis(currentObj)
        handles.KrigingAnalysisObj.plotScreeningAnalysisKrigingInterpolation(currentObj)
%         disp('OK')
    case 'Dialog Window was just closed'
%         disp('Not OK')
    otherwise
        error('Unknown Output')
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
