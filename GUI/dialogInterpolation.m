function varargout = dialogInterpolation(varargin)
% DIALOGInterpolation MATLAB code for dialogInterpolation.fig
%      DIALOGInterpolation, by itself, creates a new DIALOGInterpolation or raises the existing
%      singleton*.
%
%      H = DIALOGInterpolation returns the handle to a new DIALOGInterpolation or the handle to
%      the existing singleton*.
%
%      DIALOGInterpolation('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DIALOGInterpolation.M with the given input arguments.
%
%      DIALOGInterpolation('Property','Value',...) creates a new DIALOGInterpolation or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dialogInterpolation_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dialogInterpolation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

% Edit the above text to modify the response to help dialogInterpolation

% Last Modified by GUIDE v2.5 03-Jul-2015 13:46:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dialogInterpolation_OpeningFcn, ...
                   'gui_OutputFcn',  @dialogInterpolation_OutputFcn, ...
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


% --- Executes just before dialogInterpolation is made visible.
function dialogInterpolation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dialogInterpolation (see VARARGIN)

% Interpolation Type
strNames = varargin{2};
if varargin{1}.KrigingAnalysisObj.getnKrigingObjects==0
    msgbox('No Kriging object was selected')
    error('No Kriging object was selected')
end
handles.InterpolationType = varargin{2};

% Adjust the layout for the particular interpolation types
switch handles.InterpolationType
    case 2
        set(handles.textInVar2,'Visible','off');
        set(handles.popupmenuInVar2,'Visible','off');
        set(handles.textInVar3,'Visible','off');
        set(handles.popupmenuInVar3,'Visible','off');
%         set(handles.textInVar4,'Visible','off');
%         set(handles.popupmenuInVar4,'Visible','off');

        % Display Everything what is only used for ND-Interpolation
        set(handles.checkboxContourPlot,'Visible','off');
        set(handles.textMaxColumns,'Visible','off');
        set(handles.editNumberOfColumns,'Visible','off');
        set(handles.textCurrentRow,'Visible','off');
        set(handles.popupmenuCurrentRow,'Visible','off');
        set(handles.textMaxNumberOfRows,'Visible','off');
        set(handles.editNumberOfRows,'Visible','off');
    case 3
        set(handles.textInVar3,'Visible','off');
        set(handles.popupmenuInVar3,'Visible','off');
        
        % Display Everything what is only used for ND-Interpolation
        set(handles.checkboxContourPlot,'Visible','off');
        set(handles.textMaxColumns,'Visible','off');
        set(handles.editNumberOfColumns,'Visible','off');
        set(handles.textCurrentRow,'Visible','off');
        set(handles.popupmenuCurrentRow,'Visible','off');
        set(handles.textMaxNumberOfRows,'Visible','off');
        set(handles.editNumberOfRows,'Visible','off');
        
    case 4
        set(handles.checkboxBestData,'Visible','off');
        set(handles.checkboxConfiTube,'Visible','off');
        set(handles.checkboxContourPlot,'Visible','off');
        set(handles.checkboxDisplayOptimum,'Visible','off');
    otherwise
end

% Initialize global important variables
handles.output = hObject;
handles.KrigingAnalysisObj = varargin{1}.KrigingAnalysisObj;
handles.currentObj = get(varargin{1}.KrigingObjectivePopUp,'Value');

% Check if number of input variables is sufficient
strNames = handles.KrigingAnalysisObj.getInputVarNames(handles.currentObj);
nInputVar = handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.getnInputVar;
if nInputVar<(handles.InterpolationType-1)
    msgbox('Too less input variable for this kind of interpolation')
    error('Too less input variable for this kind of interpolation')
end

% Initialize the popup menus for the input variables
set(handles.popupmenuInVar1,'String',strNames);
set(handles.popupmenuInVar2,'String',strNames);
set(handles.popupmenuInVar3,'String',strNames);

% This has to be done even if someting different than nD_Interpolation is
% done
handles.nRows = 1;
[handles]=initializeND_Interpolation_Dialog(handles);


% Set Entries of the pop up menues
changeEntriesInRemainingPopUpMenu(handles,1);
changeEntriesInRemainingPopUpMenu(handles,2);
changeEntriesInRemainingPopUpMenu(handles,3);
changeEntriesInRemainingPopUpMenu(handles,4);

% Initialize the pop-up menu for deciding if the optimum is a minimmum or a
% maximum
contents = cellstr(get(handles.popupmenuMinOrMax,'String'));
switch handles.KrigingAnalysisObj.getMinMax(handles.currentObj)
    case -1
        newValue = find(strcmp('Minimum',contents));
    case 1
        newValue = find(strcmp('Maximum',contents));
    otherwise
        error('Unexpected value of MinMax')
end
set(handles.popupmenuMinOrMax,'Value',newValue);

% Set values of checkboxes
set(handles.checkboxConfiTube,'Value',handles.KrigingAnalysisObj.getShowBounds);
set(handles.checkboxDataPoints,'Value',handles.KrigingAnalysisObj.getShowData);
set(handles.checkboxUseInv,'Value',handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.getUseInverse);

% Set the entrie of the table
handles = setInputParametersRange(handles);
handles = setInputParametersEntries(handles);

% Update handles structure
guidata(hObject, handles);


% UIWAIT makes dialogInterpolation wait for user response (see UIRESUME)
uiwait(handles.figure1);

function [handles]=initializeND_Interpolation_Dialog(handles)
%     strNames = handles.KrigingAnalysisObj.getInputVarNames(handles.currentObj);
    nInputVar = handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.getnInputVar; 
    
    % Intialization of the indice matrix first column and second column
    % represent the x and y axis. The third column represents the index of
    % the input variable which discrete varied. The remainining columns
    % represent input variables which stays constant of the interpolation
    handles.InputVarIndiceMatrix = ones(handles.nRows,nInputVar);
    for iVar=2:nInputVar
        handles.InputVarIndiceMatrix(:,iVar)=iVar;
    end
    
    % InputVarValueMatrix saves concrete values for the input variable
    % which are used when they are fix over the interpolation. The
    % structure is similar to "InputVarIndiceMatrix"
    handles.InputVarValueMatrix = zeros(handles.nRows,nInputVar);
    
    % Set current choice of number of rows and columns in the nD plot
    set(handles.editNumberOfColumns,'String',handles.KrigingAnalysisObj.getnPlots);
    set(handles.editNumberOfRows,'String',handles.nRows);
    
    % adjust the entries in the popup for deciding which row is considered
    % right now
    changeEntriesInPopUpMenuCurrentRow(handles)
    
    % Save the names of the input variables which can be chosen be the in
    % the particular pop-up menu. Dimension is 3 X number of rows in the
    % nD-plot
    handles.InputVarStringCell = cell(3,handles.nRows);
    % Save the position of the pop-up menu. Dimension is 3 X number of rows
    % in the nD-plot
    handles.InputVarPopupValueMatrix = zeros(3,handles.nRows);
    
    

function []=changeEntriesInPopUpMenuCurrentRow(handles)

strIndicesPlots = cell(handles.nRows,1);
for iRow = 1:handles.nRows
    strIndicesPlots{iRow} = num2str(iRow);
end

set(handles.popupmenuCurrentRow,'String',strIndicesPlots);

% --- Outputs from this function are returned to the command line.
function varargout = dialogInterpolation_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;
% Get default command line output from handles structure
varargout{1} = hObject;
% varargout{2} = get(handles.popUpCoVarModel,'Value');

% The figure can be deleted now
% delete(handles.figure1);

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

function indicesWhichArechosen = getIndicesOfChosenInputVar(handles)
% Initialization
% Names of input variables
strNames = handles.KrigingAnalysisObj.getInputVarNames(handles.currentObj);

% Current chosen entry
val1 = get(handles.popupmenuInVar1,'Value');
str1 = get(handles.popupmenuInVar1,'String');

% Index of chosen input variable
index1 = strcmp(strNames,str1{val1});
index1 = find(index1==1);

% Calculate all neccessary information
indicesWhichCouldBechosen = 1:handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.getnInputVar;
indicesNotChosen = setdiff(getIndicesOfPopupInvar(handles,1),index1);
indicesWhichArechosen = setdiff(indicesWhichCouldBechosen,indicesNotChosen);

function handles = setInputParametersEntries(handles)
% This function fills the table for setting fix values of input variables
% which are not changed of the interpolation

% Initialization
strNames = handles.KrigingAnalysisObj.getInputVarNames(handles.currentObj);

% First chosen input variable
val1 = get(handles.popupmenuInVar1,'Value');
str1 = get(handles.popupmenuInVar1,'String');

% Index of chosen input variable
index1 = strcmp(strNames,str1{val1});
index1 = find(index1==1);

% Potential index of input variables
indicesWhichCouldBechosen = 1:handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.getnInputVar;

% Takes all the indices which are not yet already chosen by the popMenue
switch handles.InterpolationType
    case 2
        % Only the index in pop-up menu 1 is chosen
        indicesNotChosen = setdiff(indicesWhichCouldBechosen,val1);
    case 3
        % All indices that remain for popup Menu 1 minus the ones which is
        % actually chosen
        indicesNotChosen = setdiff(getIndicesOfPopupInvar(handles,1),index1);
    case 4
        indicesNotChosen = setdiff(getIndicesOfPopupInvar(handles,1),index1);
    otherwise
        warning('setInputParametersEntries is not defined for InterpolationType:%i',handles.InterpolationType)
end

% Input data where data points exist
dataValue = handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.getInputData;
minValue = min(dataValue);
maxValue = max(dataValue);
dataIni = cell(length(strNames)-2,2);

% Current row (only needed for nD-plot)
currentRow = get(handles.popupmenuCurrentRow,'Value');

% Collect information for the entries in the table
for iVar=1:length(strNames)-(handles.InterpolationType-1)
    dataIni{iVar,1} = strNames{indicesNotChosen(iVar)};
    
    if handles.InterpolationType==4 % nD-plot
        dataIni{iVar,2} = num2str(handles.InputVarValueMatrix(currentRow,indicesNotChosen(iVar)));
    else
        dataIni{iVar,2} = num2str(minValue(indicesNotChosen(iVar)) + 0.5*(maxValue(indicesNotChosen(iVar))-minValue(indicesNotChosen(iVar))));
    end

end

% Write information in table
set(handles.uitableInputParameters,'Data',dataIni);


function [indices]=getIndicesOfPopupInvar(handles,popMenuIndex)

% Initialization
    % Names of the input variables
strNames = handles.KrigingAnalysisObj.getInputVarNames(handles.currentObj);
    % Entries in the existing pop-up menu
val1 = get(handles.popupmenuInVar1,'Value');
str1 = get(handles.popupmenuInVar1,'String');
val2 = get(handles.popupmenuInVar2,'Value');
str2 = get(handles.popupmenuInVar2,'String');
val3 = get(handles.popupmenuInVar3,'Value');
str3 = get(handles.popupmenuInVar3,'String');
val = [val1,val2,val3];
str = {str1,str2,str3};
    % Save indices of chosen input variables
index = false(length(strNames),1);
index2 = zeros(length(val)-1,1);

runIndex = 1;
% handles.InterpolationType represents the kind of interpolation, e.g. for
% 2D interpolation only 1 variable is variated
for iIndex=setdiff(1:handles.InterpolationType-1,popMenuIndex)
    if isempty(str{iIndex})
        return
    else
        index = index|find(strcmp(str{iIndex}{val(iIndex)},strNames)==1);
    end
    % Collect the indices which are chosen in the paricular pop-up menu
    index2(runIndex) = find(strcmp(str{iIndex}{val(iIndex)},strNames)==1);
    runIndex = runIndex+1;
end
index2 = unique(index2);
indices = setdiff(1:length(strNames),index2);


% --- Executes on button press in checkboxBestData.
function checkboxBestData_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxBestData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxBestData

function []=changeEntriesInRemainingPopUpMenu(handles,ownIndex)

% Initialization
    % Names of input variables
strNames = handles.KrigingAnalysisObj.getInputVarNames(handles.currentObj);
    % handles.InterpolationType represents the kind of interpolation, e.g. for
    % 2D interpolation only 1 variable is variated
totalIndices = 1:handles.InterpolationType-1;
    % Collect the pop-up menus
popUps = {handles.popupmenuInVar1,handles.popupmenuInVar2,handles.popupmenuInVar3};

% Do work
for iIndex = totalIndices
    % Back up
    currentStr = get(popUps{iIndex},'String');
    currentStr = currentStr(get(popUps{iIndex},'Value'));
    
    % change the entries
    indices = getIndicesOfPopupInvar(handles,iIndex);
    set(popUps{iIndex},'String',strNames(indices));
    
    % only relevant for nD-plot: Remember setting when switching between
    % rows
    if sum(strcmp(currentStr,strNames(indices)) )>0
        set(popUps{iIndex},'Value',find(strcmp(currentStr,strNames(indices)) ));
    end
end
disp('')


% --- Executes on selection change in popupmenuInVar1.
function popupmenuInVar1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuInVar1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuInVar1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuInVar1

% set(handles.popupmenuInVar2,'String',{strNames{indices}});
changeEntriesInRemainingPopUpMenu(handles,1);
handles = setInputParametersEntries(handles);
handles = setInputParametersRange(handles);

if handles.InterpolationType==4
    [handles.InputVarIndiceMatrix,handles.InputVarStringCell,handles.InputVarPopupValueMatrix] = changeEntryInInputMatrices(handles,hObject,1);
end

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenuInVar1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuInVar1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuInVar2.
function popupmenuInVar2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuInVar2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuInVar2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuInVar2
changeEntriesInRemainingPopUpMenu(handles,2);
handles = setInputParametersEntries(handles);
handles = setInputParametersRange(handles);

if handles.InterpolationType==4
    [handles.InputVarIndiceMatrix,handles.InputVarStringCell,handles.InputVarPopupValueMatrix] = changeEntryInInputMatrices(handles,hObject,2);
end

% Update handles structure
guidata(hObject, handles);


function [InputVarIndiceMatrix,InputVarStringCell,InputVarPopupValueMatrix]= changeEntryInInputMatrices(handles,hObject,InVar)

% Initialiaztion
    % Number of input variables
nInputVar = handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.getnInputVar;
    % number of input variable which can be changed over the interpolation
totalIndices = 1:handles.InterpolationType-1;
    % current row of the  nD-plot
currentRow = get(handles.popupmenuCurrentRow,'Value');
    % Copy global variable for change
InputVarIndiceMatrix = handles.InputVarIndiceMatrix;
InputVarStringCell = handles.InputVarStringCell;
InputVarPopupValueMatrix = handles.InputVarPopupValueMatrix;
    % Collect existing relevant pop-up menus
popUps = {handles.popupmenuInVar1,handles.popupmenuInVar2,handles.popupmenuInVar3};

% Do work
for iIndex = totalIndices
    % All possible input variable indices which can be chosen in the
    % current selected pop-up menu
    possibleIndices = getIndicesOfPopupInvar(handles,iIndex);
    % Save current information of the chosen entry in selected pop-up menu
    InputVarIndiceMatrix(currentRow,iIndex) = ...
                                possibleIndices(get(popUps{iIndex},'Value'));
    InputVarStringCell{iIndex,currentRow} = get(popUps{iIndex},'String');
    InputVarPopupValueMatrix(iIndex,currentRow) = get(popUps{iIndex},'Value');
end
% Fill up the remaining part
InputVarIndiceMatrix(currentRow,4:end) =...
                setdiff(1:nInputVar,InputVarIndiceMatrix(currentRow,1:3));
    
    

% --- Executes during object creation, after setting all properties.
function popupmenuInVar2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuInVar2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function []=displayInputVariableBestCombination(handles,indicesRemainingInputVar)
    fprintf('The best plot(s) are create setting the non-chosen variables following values (Variable index/indices: ')
    for iVarNotChosen =indicesRemainingInputVar
        fprintf('%i\t',iVarNotChosen)
    end
    fprintf('):\n')
    handles.KrigingAnalysisObj.getChosenCombinationsForPlot


% --- Executes when entered data in editable cell(s) in uitableInputParameters.
function uitableInputParameters_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitableInputParameters (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

nInputVar = handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.getnInputVar;
currentRow = get(handles.popupmenuCurrentRow,'Value');
data = get(handles.uitableInputParameters,'Data');
counter = 1;
% for iVar=setdiff(1:nInputVar,handles.InputVarIndiceMatrix(currentRow,1:3))
if handles.InterpolationType<=3
%     for iVar=setdiff(1:nInputVar,handles.InputVarIndiceMatrix(currentRow,1:end))
%         handles.InputVarValueMatrix(currentRow,iVar) = str2double(data{counter,2});
%         counter = counter + 1;
%     end
    handles.InputVarValueMatrix(currentRow,str2double(data{counter,1})) = str2double(data{counter,2});
else
    for iVar=setdiff(1:nInputVar,handles.InputVarIndiceMatrix(currentRow,1:3))
        handles.InputVarValueMatrix(currentRow,iVar) = str2double(data{counter,2});
        counter = counter + 1;
    end
end
% handles.InputVarValueMatrix;
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbuttonInterpolation.
function pushbuttonInterpolation_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonInterpolation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nInputVar = handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.getnInputVar;
valuesRemainingInputVar = zeros(1,nInputVar-2);

val1 = get(handles.popupmenuInVar1,'Value');
val2 = get(handles.popupmenuInVar2,'Value');
% val3 = get(handles.popupmenuInVar3,'Value');
% val4 = get(handles.popupmenuInVar4,'Value');
indicesWhichCouldBechosen1 = getIndicesOfPopupInvar(handles,1);
indicesWhichCouldBechosen2 = getIndicesOfPopupInvar(handles,2);
% indicesWhichCouldBechosen3 = getIndicesOfPopupInvar(handles,3);
% indicesWhichCouldBechosen4 = getIndicesOfPopupInvar(handles,4);

% Check if the investigate range of input variables make sense
dataLB = handles.KrigingAnalysisObj.getLBInputVarInterpolation;
dataUB = handles.KrigingAnalysisObj.getUBInputVarInterpolation;
if sum(dataLB{handles.currentObj}>=dataUB{handles.currentObj})>=1
    msgbox('One of the lower bounds is bigger or equal to the associated upper bound. Please correct this before running the interpolation','Warning!')
    return
end
if get(handles.checkboxContourPlot,'Value')&&get(handles.checkboxDisplayOptimum,'Value')
    errordlg(horzcat('Uncheck either "',get(handles.checkboxContourPlot,'String'),'" or "',get(handles.checkboxDisplayOptimum,'String'),'"'))
    error(horzcat('Uncheck either "',get(handles.checkboxContourPlot,'String'),'" or "',get(handles.checkboxDisplayOptimum,'String'),'"'))
end
    
switch handles.InterpolationType
    case 2
        chosenIndices = indicesWhichCouldBechosen1(val1);

        indicesRemainingInputVar = setdiff(1:nInputVar,chosenIndices);
        data = get(handles.uitableInputParameters,'Data');
        for iVar=1:nInputVar-1
            valuesRemainingInputVar(iVar) = str2double(data{iVar,2});
        end

        
            % Do Plotting
            if get(handles.checkboxBestData,'Value')
                if get(handles.checkboxDisplayOptimum,'Value')
                    errordlg(horzcat('BestData option is not avaiable for ',get(handles.checkboxDisplayOptimum,'String')))
                else
                    handles.KrigingAnalysisObj.calcAndPlotInterpolation_2D_BestChoice(handles.currentObj,indicesWhichCouldBechosen1(val1))
                    displayInputVariableBestCombination(handles,indicesRemainingInputVar)
                end
            else
                handles.KrigingAnalysisObj.calcInterpolation_2D(handles.currentObj,indicesWhichCouldBechosen1(val1),indicesRemainingInputVar,valuesRemainingInputVar)
                if get(handles.checkboxDisplayOptimum,'Value')
                    handles.KrigingAnalysisObj.plotOptimum2D(handles.currentObj);
                else
                    handles.KrigingAnalysisObj.plotInterpolation_2D(handles.currentObj);
                end
            end
    case 3
        chosenIndices = [indicesWhichCouldBechosen1(val1),indicesWhichCouldBechosen2(val2)];

        indicesRemainingInputVar = setdiff(1:nInputVar,chosenIndices);
        data = get(handles.uitableInputParameters,'Data');
        for iVar=1:nInputVar-2
            valuesRemainingInputVar(iVar) = str2double(data{iVar,2});
        end


        % Do Plotting
        if get(handles.checkboxContourPlot,'Value')
            if get(handles.checkboxBestData,'Value')
                errordlg(horzcat('BestData option is not avaiable for ',get(handles.checkboxContourPlot,'String')))
            else
                % For Using contour plot nPlots has to be set to 1 (all other input variables are set to a defined value)
                nPlotsBackup = handles.KrigingAnalysisObj.getnPlots;
                handles.KrigingAnalysisObj.setnPlots(1);
                handles.KrigingAnalysisObj.calcInterpolation_nD(handles.currentObj,chosenIndices,indicesRemainingInputVar,valuesRemainingInputVar)
                handles.KrigingAnalysisObj.plotInterpolation_nD(handles.currentObj)
                handles.KrigingAnalysisObj.setnPlots(nPlotsBackup);
            end
        else
            if get(handles.checkboxBestData,'Value')
                if get(handles.checkboxDisplayOptimum,'Value')
                    errordlg(horzcat('BestData option is not avaiable for ',get(handles.checkboxDisplayOptimum,'String')))
                else
                    handles.KrigingAnalysisObj.calcAndPlotInterpolation_3D_BestChoice(handles.currentObj,chosenIndices)
                    displayInputVariableBestCombination(handles,indicesRemainingInputVar)
                end
            else
                handles.KrigingAnalysisObj.calcInterpolation_3D(handles.currentObj,chosenIndices,indicesRemainingInputVar,valuesRemainingInputVar)
                if get(handles.checkboxDisplayOptimum,'Value')
                    handles.KrigingAnalysisObj.plotOptimum3D(handles.currentObj);
                else
                    handles.KrigingAnalysisObj.plotInterpolation_3D(handles.currentObj);
                end
            end
        end
        
        
    case 4
       
        valuesRemainingInputVarValues = zeros(handles.nRows,nInputVar-3);
        for iVar=1:handles.nRows
            valuesRemainingInputVarValues(iVar,:) = handles.InputVarValueMatrix(iVar,handles.InputVarIndiceMatrix(iVar,4:end));
        end

        % Do Plotting
        handles.KrigingAnalysisObj.calcInterpolation_nD(handles.currentObj,handles.InputVarIndiceMatrix(:,1:3),handles.InputVarIndiceMatrix(:,4:end),valuesRemainingInputVarValues);
        handles.KrigingAnalysisObj.plotInterpolation_nD(handles.currentObj);

    otherwise
        msgbox('Interpolation Type is not defined')
end



% --- Executes on selection change in popupmenuInVar3.
function popupmenuInVar3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuInVar3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuInVar3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuInVar3
changeEntriesInRemainingPopUpMenu(handles,3);
handles = setInputParametersEntries(handles);
handles = setInputParametersRange(handles);

if handles.InterpolationType==4
    [handles.InputVarIndiceMatrix,handles.InputVarStringCell,handles.InputVarPopupValueMatrix] = changeEntryInInputMatrices(handles,hObject,3);
end

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenuInVar3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuInVar3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuInVar4.
function popupmenuInVar4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuInVar4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuInVar4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuInVar4
changeEntriesInRemainingPopUpMenu(handles,4);
handles = setInputParametersEntries(handles);
handles = setInputParametersRange(handles);

% --- Executes during object creation, after setting all properties.
function popupmenuInVar4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuInVar4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Handle of manually created figures have integers but the GUI figure have
% float-values handles
fh=findall(0,'type','figure');
ver = version('-release');
if str2double(ver(1:4))>2014
    indices = [fh.Number]';
else
    indices = fh;
end
close(indices((indices-round(indices))==0))



% --- Executes on button press in checkboxConfiTube.
function checkboxConfiTube_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxConfiTube (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% set(handles.checkboxConfiTube,'Value',handles.KrigingAnalysisObj.getShowBounds)
% set(handles.checkboxDataPoints,'Value',handles.KrigingAnalysisObj.getShowData)
% Hint: get(hObject,'Value') returns toggle state of checkboxConfiTube
handles.KrigingAnalysisObj.setShowBounds(logical(get(hObject,'Value')))
guidata(hObject, handles);

% --- Executes on button press in checkboxDataPoints.
function checkboxDataPoints_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxDataPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkboxDataPoints
handles.KrigingAnalysisObj.setShowData(logical(get(hObject,'Value')))
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function textInVar1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textInVar1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --------------------------------------------------------------------
function uitableMinMaxVar_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to uitableMinMaxVar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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

indicesOfChosenInputVar = getIndicesOfChosenInputVar(handles);

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
    
%     handles
end


% --- Executes on button press in checkboxContourPlot.
function checkboxContourPlot_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxContourPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxContourPlot


% --- Executes on button press in checkboxUseInv.
function checkboxUseInv_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxUseInv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxUseInv
% set(handles.checkboxUseInv,'Value',handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.getUseInverse);
handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.setUseInverse(logical(get(hObject,'Value')))


% --- Executes on button press in checkboxDisplayOptimum.
function checkboxDisplayOptimum_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxDisplayOptimum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxDisplayOptimum


% --- Executes on selection change in popupmenuMinOrMax.
function popupmenuMinOrMax_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuMinOrMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuMinOrMax contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuMinOrMax
contents = cellstr(get(hObject,'String'));
switch contents{get(hObject,'Value')}
    case 'Minimum'
        handles.KrigingAnalysisObj.setMinMax(handles.currentObj,-1)
    case 'Maximum'
        handles.KrigingAnalysisObj.setMinMax(handles.currentObj,1)
    otherwise
        errordlg('Unexpected Choice')
end
    
    

% --- Executes during object creation, after setting all properties.
function popupmenuMinOrMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuMinOrMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editNumberOfColumns_Callback(hObject, eventdata, handles)
% hObject    handle to editNumberOfColumns (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNumberOfColumns as text
%        str2double(get(hObject,'String')) returns contents of editNumberOfColumns as a double
handles.KrigingAnalysisObj.setnPlots(floor(str2double(get(hObject,'String'))))
set(hObject,'String',num2str(handles.KrigingAnalysisObj.getnPlots))
changeEntriesInPopUpMenuCurrentRow(handles)



% --- Executes during object creation, after setting all properties.
function editNumberOfColumns_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNumberOfColumns (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuInVar4.
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuInVar4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuInVar4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuInVar4


% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuInVar4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function textCurrentRow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textCurrentRow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in popupmenuCurrentRow.
function popupmenuCurrentRow_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuCurrentRow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuCurrentRow contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuCurrentRow

% Initialization
totalIndices = 1:handles.InterpolationType-1;
currentRow = get(handles.popupmenuCurrentRow,'Value');
popUps = {handles.popupmenuInVar1,handles.popupmenuInVar2,handles.popupmenuInVar3};

for iIndex = totalIndices
    if ~isempty(handles.InputVarStringCell{iIndex,currentRow})
        set(popUps{iIndex},'String',handles.InputVarStringCell{iIndex,currentRow})
        set(popUps{iIndex},'Value',handles.InputVarPopupValueMatrix(iIndex,currentRow))
    end
end

handles = setInputParametersEntries(handles);
handles = setInputParametersRange(handles);


guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function popupmenuCurrentRow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuCurrentRow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function uipanel1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function editNumberOfRows_Callback(hObject, eventdata, handles)
% hObject    handle to editNumberOfRows (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNumberOfRows as text
%        str2double(get(hObject,'String')) returns contents of editNumberOfRows as a double
handles.nRows = floor(str2double(get(hObject,'String')));
set(hObject,'String',num2str(handles.nRows))

[handles]=initializeND_Interpolation_Dialog(handles);
changeEntriesInPopUpMenuCurrentRow(handles)
handles = setInputParametersEntries(handles);
handles = setInputParametersRange(handles); 

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function editNumberOfRows_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNumberOfRows (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function uitableMinMaxVar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitableMinMaxVar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
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
