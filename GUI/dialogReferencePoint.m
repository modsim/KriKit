function varargout = dialogReferencePoint(varargin)
% DIALOGREFERENCEPOINT MATLAB code for dialogReferencePoint.fig
%      DIALOGREFERENCEPOINT, by itself, creates a new DIALOGREFERENCEPOINT or raises the existing
%      singleton*.
%
%      H = DIALOGREFERENCEPOINT returns the handle to a new DIALOGREFERENCEPOINT or the handle to
%      the existing singleton*.
%
%      DIALOGREFERENCEPOINT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DIALOGREFERENCEPOINT.M with the given input arguments.
%
%      DIALOGREFERENCEPOINT('Property','Value',...) creates a new DIALOGREFERENCEPOINT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dialogReferencePoint_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dialogReferencePoint_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% Edit the above text to modify the response to help dialogReferencePoint
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

% Last Modified by GUIDE v2.5 03-Jul-2015 13:48:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dialogReferencePoint_OpeningFcn, ...
                   'gui_OutputFcn',  @dialogReferencePoint_OutputFcn, ...
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


% --- Executes just before dialogReferencePoint is made visible.
function dialogReferencePoint_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dialogReferencePoint (see VARARGIN)

if varargin{1}.KrigingAnalysisObj.getnKrigingObjects==0
    msgbox('No Kriging object was selected')
    error('No Kriging object was selected')
end

% Choose default command line output for dialogCovarEsti
handles.output = hObject;
handles.KrigingAnalysisObj = varargin{1}.KrigingAnalysisObj;
handles.currentObj = get(varargin{1}.KrigingObjectivePopUp,'Value');

% Set the entrie of the table
handles = setInputParametersRange(handles);
handles=setSovlerTableEntries(handles);

% Choose default command line output for dialogReferencePoint
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dialogReferencePoint wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = dialogReferencePoint_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;
% Get default command line output from handles structure
% varargout{1} = handles.KrigingAnalysisObj;
if isfield(handles,'pushbutton1')
    varargout{1} = 'OK Button Was Pushed';
    delete(handles.figure1);
else
    varargout{1} = 'Dialog Window was just closed';
end
% varargout{2} = get(handles.popUpCoVarModel,'Value');

% The figure can be deleted now

% close(hObject);

% % --- Executes when user attempts to close figure1.
% function figure1_CloseRequestFcn(hObject, eventdata, handles)
% % hObject    handle to figure1 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hint: delete(hObject) closes the figure
% if isequal(get(hObject,'waitstatus'),'waiting')
%     uiresume(hObject);
% else
%     delete(hObject);
% end

function handles = setSovlerTableEntries(handles)
nInputVar = handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.getnInputVar;
inputVarNames = handles.KrigingAnalysisObj.getInputVarNames(handles.currentObj);
refPoints = handles.KrigingAnalysisObj.getReferencePoint;
if isempty(refPoints)||size(refPoints,2)~=nInputVar
    refPoints=zeros(1,nInputVar);
    handles.KrigingAnalysisObj.setReferencePoint(refPoints);
end
data=cell(nInputVar,2);

for iInputVar = 1:nInputVar
    data{iInputVar,1}=inputVarNames{iInputVar};
    data{iInputVar,2}=num2str(refPoints(iInputVar));  
end
set(handles.uitableReferenceValues,'Data',data)

% --- Executes when entered data in editable cell(s) in uitableReferenceValues.
function uitableReferenceValues_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitableReferenceValues (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
% data = get(handles.uitableSolverSettings,'Data');
newValue = str2double(eventdata.NewData);
refPoints = handles.KrigingAnalysisObj.getReferencePoint;
refPoints(eventdata.Indices(1)) = newValue;
handles.KrigingAnalysisObj.setReferencePoint(refPoints);

% % --- Executes during object creation, after setting all properties.
% function figure1_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to figure1 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% 
% % --- Executes during object deletion, before destroying properties.
% function figure1_DeleteFcn(hObject, eventdata, handles)
% % hObject    handle to figure1 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)


function handles = setInputParametersRange(handles)

% Initialization
    % Names of input variables
strNames = handles.KrigingAnalysisObj.getInputVarNames(handles.currentObj);
    % Allowed ranges of input variable values for the interpolation 
dataLBRange = handles.KrigingAnalysisObj.getLBInputVarInterpolation();
dataUBRange = handles.KrigingAnalysisObj.getUBInputVarInterpolation();
dataLBRange = dataLBRange{handles.currentObj};
dataUBRange = dataUBRange{handles.currentObj};
dataIni = cell(length(strNames),2);
    % If no range is adjusted yet, take min and max of the provided data
if isempty(dataLBRange)||isempty(dataUBRange)
    dataLBRange = min(handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.getInputData);
    dataUBRange = max(handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.getInputData);
    % Save for the future
    handles.KrigingAnalysisObj.setLBInputVarInterpolation(handles.currentObj,dataLBRange);
    handles.KrigingAnalysisObj.setUBInputVarInterpolation(handles.currentObj,dataUBRange);
end
    % Write in table
for iVar=1:length(strNames)
    iVarIndex = iVar;
    dataIni{iVar,1} = strNames{iVarIndex};
    dataIni{iVar,2} = num2str(dataLBRange(iVarIndex));
    dataIni{iVar,3} = num2str(dataUBRange(iVarIndex));
end
    % Save
set(handles.uitableMinMaxVar,'Data',dataIni);

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isequal(get(handles.figure1,'waitstatus'),'waiting')
    uiresume(handles.figure1);
else
    delete(handles.figure1);
end


% --- Executes when entered data in editable cell(s) in uitableMinMaxVar.
function uitableMinMaxVar_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitableMinMaxVar (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
columns = get(hObject,'ColumnName');
data = get(hObject,'Data');

newValue = str2double(eventdata.NewData);
if ~isnan(newValue)
    switch columns{eventdata.Indices(2)}
        case 'Min'
            handles.KrigingAnalysisObj.setLBInputVarInterpolation(handles.currentObj,str2double(data(:,2)))
        case 'Max'
            handles.KrigingAnalysisObj.setUBInputVarInterpolation(handles.currentObj,str2double(data(:,3)))
        otherwise
            errror('Unexpectd column name of the edited column in uitableCovarPara')
    end
    
    if sum(str2double(data(:,2))>=str2double(data(:,3)))
        msgbox('One of the lower bounds is bigger or equal to the associated upper bound. Please correct this before running the interpolation','Warning!')
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
