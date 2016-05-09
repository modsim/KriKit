function varargout = dialogCovarEsti(varargin)
% DIALOGCOVARESTI MATLAB code for dialogCovarEsti.fig
%      DIALOGCOVARESTI, by itself, creates a new DIALOGCOVARESTI or raises the existing
%      singleton*.
%
%      H = DIALOGCOVARESTI returns the handle to a new DIALOGCOVARESTI or the handle to
%      the existing singleton*.
%
%      DIALOGCOVARESTI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DIALOGCOVARESTI.M with the given input arguments.
%
%      DIALOGCOVARESTI('Property','Value',...) creates a new DIALOGCOVARESTI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dialogCovarEsti_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dialogCovarEsti_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% Edit the above text to modify the response to help dialogCovarEsti
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

% Last Modified by GUIDE v2.5 11-Sep-2015 15:11:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dialogCovarEsti_OpeningFcn, ...
                   'gui_OutputFcn',  @dialogCovarEsti_OutputFcn, ...
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


% --- Executes just before dialogCovarEsti is made visible.
function dialogCovarEsti_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dialogCovarEsti (see VARARGIN)

if varargin{1}.KrigingAnalysisObj.getnKrigingObjects==0
    msgbox('No Kriging object was selected')
    error('No Kriging object was selected')
end

% Choose default command line output for dialogCovarEsti
handles.output = hObject;
handles.KrigingAnalysisObj = varargin{1}.KrigingAnalysisObj;
handles.currentObj = get(varargin{1}.KrigingObjectivePopUp,'Value');

% Set popup menues to current Kriging setting
set(handles.popUpVargramEstiType,'Value',handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.getCovariogramEstimationType)
set(handles.popUpCoVarModel,'Value',handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.getCovariogramModelChoice)
set(handles.popUpChooseSolver,'Value',handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.getUseSolver)

% Set Text Box
handles = setFormulaBox(handles,true);

% Set the entries in the table containing the solver settings 
handles = setSovlerTableEntries(handles);
handles = setCoVarParaTableEntries(handles);

% check if estimation is still running
handles.checkIfEstimationIsRunning = 0;

% Set Checkbox
set(handles.checkboxDetails,'Value',logical(handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.getShowDetails()));


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dialogCovarEsti wait for user response (see UIRESUME)
uiwait(handles.figure1);

function handles = setFormulaBox(handles,firstTimeBool)

switch handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.getCovariogramModelChoice
    case 1
        c1 = '$$\mathrm{cov}\left( \mathbf{x}_{i},\mathbf{x}_{j} \right) = \sigma^2\exp\left( - \left( \frac{0.5r}{2\theta^2} \right)\right)$$';
        c2 = '$$r=\left( \mathbf{x}_{i} - \mathbf{x}_{j} \right)^{T}\left( \mathbf{x}_{i} - \mathbf{x}_{j} \right)$$';
        chosenFormula = {c1,'',c2};
%         chosenFormula = '$$\mathrm{cov}\left( \mathbf{x}_{i},\mathbf{x}_{j} \right) = \sigma^2\exp\left( - \left( \frac{\sqrt{\left( \mathbf{x}_{i} - \mathbf{x}_{j} \right)^{T}\left( \mathbf{x}_{i} - \mathbf{x}_{j} \right)}}{2\theta^2} \right)\right)$$';
    case 2
        c1 = '$$\mathrm{cov}\left( \mathbf{x}_{i},\mathbf{x}_{j} \right) = \sigma_{{Nugget}}^2 + \sigma^2\exp\left( - \left( \frac{0.5r}{2\theta^2} \right)\right)$$';
        c2 = '$$r=\left( \mathbf{x}_{i} - \mathbf{x}_{j} \right)^{T}\left( \mathbf{x}_{i} - \mathbf{x}_{j} \right)$$';
        chosenFormula = {c1,'',c2};
    case 3
        c1 = '$$\mathrm{cov}\left( \mathbf{x}_{i},\mathbf{x}_{j} \right) = \sigma_{{Nugget}}^2 + \sigma^2\exp\left( - \left( \frac{0.5r}{2\theta^2} \right)\right)$$';
        c2 = '$$r = \sqrt{\sum_{l = 1}^{k}\frac{\left( x_{i,l} - x_{j,l} \right)^{p_{l}}}{\theta_{l}^{2}}}$$';
        chosenFormula = {c1,'',c2};
%         chosenFormula = '$$\mathrm{cov}\left( \mathbf{x}_{i},\mathbf{x}_{j} \right) = \sigma_{Nugget}^2 + \sigma^2\exp\left( - \left( \sum_{l = 1}^{k}{\frac{\left( x_{i,l} - x_{j,l} \right)^{p_{l}}}{\theta_{l}^2}} \right) \right)$$';
    case 4
        c1 = '$$\mathrm{cov}\left( \mathbf{x}_{i},\mathbf{x}_{j} \right) = \sigma_{{Nugget}}^2 + \sigma^2\exp\left( - \left( \frac{0.5r}{2\theta^2} \right)\right)$$';
        c2 = '$$r = \sqrt{\sum_{l = 1}^{k}\frac{\left( x_{i,l} - x_{j,l} \right)^{2}}{\theta_{l}^{2}}}$$';
        chosenFormula = {c1,'',c2};
    case 5
        c = '$$\mathrm{cov}\left( \mathbf{x}_{i},\mathbf{x}_{j} \right) = \sigma_{{Nugget}}^{2} + \sigma^{2}\left( 1 + \frac{\left( \sqrt{3}r \right)}{\theta^{2}} \right)\exp\left( - \frac{\left( \sqrt{3}r \right)}{\theta^{2}} \right),$$';
        c2 = '$$r = \sqrt{\left( \mathbf{x}_{i} - \mathbf{x}_{j} \right)^{T}\left( \mathbf{x}_{i} - \mathbf{x}_{j} \right)}$$';
        chosenFormula = {c,c2};
    case 6
        c = '$$\mathrm{cov}\left( \mathbf{x}_{i},\mathbf{x}_{j} \right) = \sigma_{Nugget}^{2} + \sigma^{2}\left( 1 + \sqrt{3}r \right)\exp\left( - \sqrt{3}r \right)$$';
        c2 = '$$r = \sqrt{\sum_{l = 1}^{k}\frac{\left( x_{i,l} - x_{j,l} \right)^{2}}{\theta_{l}^{2}}}$$';
        chosenFormula = {c,c2};
    case 7
        c = '$$\mathrm{cov}\left( \mathbf{x}_{i},\mathbf{x}_{j} \right) = \sigma_{Nugget}^{2} + \sigma^{2}\left( 1 + \frac{\sqrt{5}r}{\theta^{2}} + \frac{5r^{2}}{3\theta^{2}} \right)\exp\left( - \frac{\sqrt{5}r}{\theta^{2}} \right),$$';
        c2 = '$$r = \ \sqrt{\left( \mathbf{x}_{i} - \mathbf{x}_{j} \right)^{T}\left( \mathbf{x}_{i} - \mathbf{x}_{j} \right)}$$';
        chosenFormula = {c,c2};
    case 8
        c = '$$\mathrm{cov}\left( \mathbf{x}_{i},\mathbf{x}_{j} \right) = \sigma_{Nugget}^{2} + \sigma^{2}\left( 1 + \sqrt{3}r + \frac{5}{3}r^{2} \right)\exp\left( - \sqrt{5}r \right)$$';
        c2 = '$$r = \sqrt{\sum_{l = 1}^{k}\frac{\left( x_{i,l} - x_{j,l} \right)^{2}}{\theta_{l}^{2}}}$$';
        chosenFormula = {c,c2};
    otherwise
        chosenFormula = 'Not defined';
end
        
hObjectHere=handles.textBoxCovModel;
% lbls = findobj(handles,'-regexp','tag','textBoxCovModel');

% Get current text, position and tag
set(hObjectHere,'units','normalized');
% stringFormula = get(hObjectHere,'string');
positionFormula = get(hObjectHere,'position');
tagSaved = get(hObjectHere,'tag');

% Remove the UICONTROL
% if firstTimeBool
delete(hObjectHere);
% end

% Replace it with a TEXT object 
handles.textBoxCovModel = text(positionFormula(1),positionFormula(2),chosenFormula,'interpreter','latex','FontSize',20);

axis off; 

function handles = setSovlerTableEntries(handles)
val = get(handles.popUpChooseSolver,'Value');
str = get(handles.popUpChooseSolver,'String');
switch str{val}
    case 'Fminsearch'
        data = {'Max. Number of Iterations',num2str(handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.getnIterationsSolver)};
    case 'Fmincon'
        data = {'Max. Number of Iterations',num2str(handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.getnIterationsSolver)};
    case 'Open Source Genetic Algorithm'
        data = {'Max. Number of Generations',num2str(handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.getGenerations);...
                'Time Limit in Sec.',num2str(handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.getTimeLimit);...
                'Number of Individuals in Population',num2str(handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.getPopulationSize)};
    case 'Optimization Toolbox Genetic Algorithm'
        data = {'Max. Number of Generations',num2str(handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.getGenerations);...
                'Time Limit in Sec.',num2str(handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.getTimeLimit);...
                'Number of Individuals in Population',num2str(handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.getPopulationSize)};
    case 'Matlabs Gaussian Process Regression'
        data = {};
    otherwise
        error('%s is not an option for the solver',str{val})
end
set(handles.uitableSolverSettings,'Data',data)


function handles = setCoVarParaTableEntries(handles)
val = get(handles.popUpCoVarModel,'Value');
str = get(handles.popUpCoVarModel,'String');
LBparaValues = handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.getLBCovariogramModelParameters;
UBparaValues = handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.getUBCovariogramModelParameters;
para = handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.getCovariogramModelParameters;
nInputVar = handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.getnInputVar;

switch str{val}
    case 'Squared Exponential'
        if ~isempty(LBparaValues)
            data = {'theta',num2str(LBparaValues(1)),num2str(UBparaValues(1)),num2str(para(1));...
                    'sigma',num2str(LBparaValues(2)),num2str(UBparaValues(2)),num2str(para(2))};
        else
            data = {'theta',num2str(-inf),num2str(inf),num2str(para(1));...
                    'sigma',num2str(-inf),num2str(inf),num2str(para(2))};
        end
    case 'Squared Exponential with Nugget'
        if ~isempty(LBparaValues)
            data = {'theta',num2str(LBparaValues(1)),num2str(UBparaValues(1)),num2str(para(1));...
                    'sigma',num2str(LBparaValues(2)),num2str(UBparaValues(2)),num2str(para(2));...
                    'nugget',num2str(LBparaValues(3)),num2str(UBparaValues(3)),num2str(para(3))};
        else
            data = {'theta',num2str(-inf),num2str(inf),num2str(para(1));...
                    'sigma',num2str(-inf),num2str(inf),num2str(para(2));...
                    'nugget',num2str(-inf),num2str(inf),num2str(para(3))};
        end
    case 'p-Exponential Function with separate length scale'
        % Special Case
        data  = cell(2*nInputVar+2,3);
        
        % First the stings
        for iVar = 1:nInputVar
            data{iVar,1} = horzcat('theta',num2str(iVar));
            data{iVar+nInputVar,1} = horzcat('p',num2str(iVar));
        end
        data{end-1,1} = 'sigma';
        data{end,1} = 'nugget';
        data{end-1,4} = num2str(para(end-1));
        data{end,4} = num2str(para(end));
        
        % Than the LB/UB-Values
        if length(LBparaValues)==2*nInputVar+2
            for iVar = 1:nInputVar
                data{iVar,2} = num2str(LBparaValues(iVar));
                data{iVar+nInputVar,2} = num2str(LBparaValues(iVar+nInputVar));
            end
            data{end-1,2} = num2str(LBparaValues(end-1));
            data{end,2} = num2str(LBparaValues(end));
        else
            for iVar = 1:nInputVar
                data{iVar,2} = num2str(-inf);
                data{iVar+nInputVar,2} = num2str(-inf);
            end
            data{end-1,2} = num2str(-inf);
            data{end,2} = num2str(-inf);
        end
        if length(UBparaValues)==2*nInputVar+2
            for iVar = 1:nInputVar
                data{iVar,3} = num2str(UBparaValues(iVar));
                data{iVar+nInputVar,3} = num2str(UBparaValues(iVar+nInputVar));
            end
            data{end-1,3} = num2str(UBparaValues(end-1));
            data{end,3} = num2str(UBparaValues(end));
        else
            for iVar = 1:nInputVar
                data{iVar,3} = num2str(inf);
                data{iVar+nInputVar,3} = num2str(inf);
            end
            data{end-1,3} = num2str(inf);
            data{end,3} = num2str(inf);
        end

        % Current parameter Values
        for iVar = 1:nInputVar
            data{iVar,4} = num2str(para(iVar));
            data{iVar+nInputVar,4} = num2str(para(iVar+nInputVar));
        end
        
    case 'Squared Exponential Kernel with separate length scale'
        data  = cell(nInputVar+2,3);
        
%         % First the stings
%         for iVar = 1:nInputVar
%             data{iVar,1} = horzcat('theta',num2str(iVar));
%         end
%         data{end-1,1} = 'sigma';
%         data{end,1} = 'nugget';
%         data{end-1,3} = para(end-1);
%         data{end,3} = para(end);
        
        [data] = setParameterBoundInTable(data,LBparaValues,UBparaValues,nInputVar,para);
        
    case 'Matern 3/2' 
        if ~isempty(LBparaValues)
            data = {'theta',num2str(LBparaValues(1)),num2str(UBparaValues(1)),num2str(para(1));...
                    'sigma',num2str(LBparaValues(2)),num2str(UBparaValues(2)),num2str(para(2));...
                    'nugget',num2str(LBparaValues(3)),num2str(UBparaValues(3)),num2str(para(3))};
        else
            data = {'theta',num2str(-inf),num2str(inf),num2str(para(1));...
                    'sigma',num2str(-inf),num2str(inf),num2str(para(2));...
                    'nugget',num2str(-inf),num2str(inf),num2str(para(3))};
        end
    case 'Matern 3/2 with separate length scale'
        data  = cell(nInputVar+2,3);
        
        % First the stings
        for iVar = 1:nInputVar
            data{iVar,1} = horzcat('theta',num2str(iVar));
        end
        data{end-1,1} = 'sigma';
        data{end,1} = 'nugget';
        
        [data] = setParameterBoundInTable(data,LBparaValues,UBparaValues,nInputVar,para);
    case 'Matern 5/2'
        if ~isempty(LBparaValues)
            data = {'theta',num2str(LBparaValues(1)),num2str(UBparaValues(1)),num2str(para(1));...
                    'sigma',num2str(LBparaValues(2)),num2str(UBparaValues(2)),num2str(para(2));...
                    'nugget',num2str(LBparaValues(3)),num2str(UBparaValues(3)),num2str(para(3))};
        else
            data = {'theta',num2str(-inf),num2str(inf),num2str(para(1));...
                    'sigma',num2str(-inf),num2str(inf),num2str(para(2));...
                    'nugget',num2str(-inf),num2str(inf),num2str(para(3))};
        end
    case 'Matern 5/2 with separate length scale'
        data  = cell(nInputVar+2,4);
        
        % First the stings
        for iVar = 1:nInputVar
            data{iVar,1} = horzcat('theta',num2str(iVar));
        end
        data{end-1,1} = 'sigma';
        data{end,1} = 'nugget';
        
        [data] = setParameterBoundInTable(data,LBparaValues,UBparaValues,nInputVar,para);
    otherwise
end
set(handles.uitableCovarPara,'Data',data)

function [data] = setParameterBoundInTable(data,LBparaValues,UBparaValues,nInputVar,para)
    % First the stings
    for iVar = 1:nInputVar
        data{iVar,1} = horzcat('theta',num2str(iVar));
    end
    data{end-1,1} = 'sigma';
    data{end,1} = 'nugget';
    data{end-1,4} = num2str(para(end-1));
    data{end,4} = num2str(para(end));

    % Than the LB/UB-Values
    if length(LBparaValues)==nInputVar+2
        for iVar = 1:nInputVar
            data{iVar,2} = num2str(LBparaValues(iVar));
%             data{iVar+nInputVar,2} = num2str(LBparaValues(iVar+nInputVar));
        end
        data{end-1,2} = num2str(LBparaValues(end-1));
        data{end,2} = num2str(LBparaValues(end));
    else
        for iVar = 1:nInputVar
            data{iVar,2} = num2str(-inf);
%             data{iVar+nInputVar,2} = num2str(-inf);
        end
        data{end-1,2} = num2str(-inf);
        data{end,2} = num2str(-inf);
    end
    if length(UBparaValues)==nInputVar+2
        for iVar = 1:nInputVar
            data{iVar,3} = num2str(UBparaValues(iVar));
%             data{iVar+nInputVar,3} = num2str(UBparaValues(iVar+nInputVar));
        end
        data{end-1,3} = num2str(UBparaValues(end-1));
        data{end,3} = num2str(UBparaValues(end));
    else
        for iVar = 1:nInputVar
            data{iVar,3} = num2str(inf);
%             data{iVar+nInputVar,3} = num2str(inf);
        end
        data{end-1,3} = num2str(inf);
        data{end,3} = num2str(inf);
    end
    
    % Current parameter Values
    for iVar = 1:nInputVar
        data{iVar,4} = num2str(para(iVar));
%         data{iVar+nInputVar,4} = num2str(para(iVar+nInputVar));
    end

    
% --- Outputs from this function are returned to the command line.
function varargout = dialogCovarEsti_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.KrigingAnalysisObj;
% varargout{2} = get(handles.popUpCoVarModel,'Value');

% The figure can be deleted now
delete(handles.figure1);


% --- Executes on selection change in popUpVargramEstiType.
function popUpVargramEstiType_Callback(hObject, eventdata, handles)
% hObject    handle to popUpVargramEstiType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popUpVargramEstiType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popUpVargramEstiType
handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.setCovariogramEstimationType(get(handles.popUpVargramEstiType,'Value'))

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popUpVargramEstiType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popUpVargramEstiType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isequal(get(hObject,'waitstatus'),'waiting')
    uiresume(hObject);
else
    delete(hObject);
end


% --- Executes on selection change in popUpCoVarModel.
function popUpCoVarModel_Callback(hObject, eventdata, handles)
% hObject    handle to popUpCoVarModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popUpCoVarModel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popUpCoVarModel

% Set new Value
handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.setCovariogramModelChoice(get(handles.popUpCoVarModel,'Value'))

% Set the entries in the table containing the solver settings 
handles = setCoVarParaTableEntries(handles);

% Reset Formula
handles = setFormulaBox(handles,false);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popUpCoVarModel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popUpCoVarModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxDetails.
function checkboxDetails_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxDetails (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxDetails
handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.setShowDetails(get(handles.checkboxDetails,'Value'));

% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in popUpChooseSolver.
function popUpChooseSolver_Callback(hObject, eventdata, handles)
% hObject    handle to popUpChooseSolver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popUpChooseSolver contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popUpChooseSolver
% Set new Value

% Initialization
stringCell = get(handles.popUpChooseSolver,'String');
% Set setUseMatlabRegressionGP by default off
handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.setUseMatlabRegressionGP(false)

switch stringCell{get(handles.popUpChooseSolver,'Value')}
    case 'Fminsearch'
        handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.setUseSolver(1)
    case 'Fmincon'
        handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.setUseSolver(2)
    case 'Open Source Genetic Algorithm'
        handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.setUseSolver(3)
    case 'Optimization Toolbox Genetic Algorithm'
        handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.setUseSolver(4)
    case 'Matlabs Gaussian Process Regression'
        necessaryFeatures = {'Statistics_Toolbox'};
        validLincence = cellfun(@(f) license('checkout',f),necessaryFeatures);
        versionStr = version;

        if ~validLincence||(str2double(versionStr(end-5:end-2))<2015)||(str2double(versionStr(end-5:end-2))==2015&&strcmp(versionStr(end-1),'a'))
            errordlg('For Using RegressionGP, you need at least Matlab version 2015b and the Statistic Toolbox')
        else
            handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.setUseMatlabRegressionGP(true)
        end
    otherwise
        errordlg('Solver unknown')
end

% Set the entries in the table containing the solver settings 
handles = setSovlerTableEntries(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popUpChooseSolver_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popUpChooseSolver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when entered data in editable cell(s) in uitableSolverSettings.
function uitableSolverSettings_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitableSolverSettings (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
data = get(handles.uitableSolverSettings,'Data');
newValue = str2double(eventdata.NewData);

if ~isnan(newValue)

    switch data{eventdata.Indices(1),1}
        case 'Max. Number of Iterations'
            handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.setnIterationsSolver(round(abs(newValue)));
        case 'Max. Number of Generations'
            handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.setGenerations(round(abs(newValue)));
%             fprintf('Gen was changed to %d \n',round(abs(newValue)))
        case 'Time Limit in Sec.'
            handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.setTimeLimit(round(abs(newValue)));
%             fprintf('Time was changed to %d \n',round(abs(newValue)))
        case 'Number of Individuals in Population'
            handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.setPopulationSize(round(abs(newValue)));
%             fprintf('Pop was changed to %d \n',round(abs(newValue)))
        otherwise
            error('%s is not in the list of allowed settings for the solver',data{eventdata.Indices(1),1})
    end
    
end

% Update handles structure
handles=setSovlerTableEntries(handles);
guidata(hObject, handles);


% --- Executes when entered data in editable cell(s) in uitableCovarPara.
function uitableCovarPara_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitableCovarPara (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
columns = get(handles.uitableCovarPara,'ColumnName');
data = get(handles.uitableCovarPara,'Data');
newValue = str2double(eventdata.NewData);
if ~isnan(newValue)
    switch columns{eventdata.Indices(2)}
        case 'Lower Bound'
            handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.setLBCovariogramModelParameters( str2double({data{:,eventdata.Indices(2)}}) )
        case 'Upper Bound'
            handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.setUBCovariogramModelParameters( str2double({data{:,eventdata.Indices(2)}}) )
        case 'Current Value'
            handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.setCovariogramModelParameters( str2double({data{:,eventdata.Indices(2)}}) )
            handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.setInitialCovariogramParameters( str2double({data{:,eventdata.Indices(2)}}) )
        otherwise
            errror('Unexpectd column name of the edited column in uitableCovarPara')
    end
%     handles
end

% Update handles structure
handles=setCoVarParaTableEntries(handles);
guidata(hObject, handles);


% --- Executes on button press in pushbuttonEstimatePara.
function pushbuttonEstimatePara_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonEstimatePara (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Make sure that drawNow has an effect (further below)
set(hObject,'Interruptible','off')

val = get(handles.popUpChooseSolver,'Value');
str = get(handles.popUpChooseSolver,'String');
if strcmp('Fmincon',str{val})||strcmp('Fminsearch',str{val});
    data = get(handles.uitableCovarPara,'Data');
    handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.setInitialCovariogramParameters( str2double({data{:,4}}) );
%     handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.setInitialCovariogramParameters(handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.getCovariogramModelParameters);
end


mBox1=msgbox('Optimization Starts, new dialog is opened as soon as optimization finished');
% mBox1=questdlg('Optimization Starts, new dialog is opened as soon as optimization finished','Optimization Runs','OK','Abort','OK2');
switch str{val}
    % Use FminSearch if no Optimization_Toolbox is available how ever this
    % algorithm does not considers constraints!
    case 'Fminsearch'
        handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.estimateVariance;
        handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.solveLeastSquareCovariogramFminSearch;
    case 'Fmincon'

        if (license('checkout','Optimization_Toolbox'))
            handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.estimateVariance;
            handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.solveLeastSquareCovariogram;
        else
            msgbox('No Optimization Toolbox licence available')
        end

    case 'Open Source Genetic Algorithm'
        handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.estimateVariance;
        handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.solveLeastSquareCovariogramGA2;
    case 'Optimization Toolbox Genetic Algorithm'

        if (license('checkout','GADS_Toolbox'))
            handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.estimateVariance;
            handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.solveLeastSquareCovariogramGA;
        else
            msgbox('No Global Optimization Toolbox licence available')
        end
    case 'Matlabs Gaussian Process Regression'
        testIndices = [2,4:8];
        modelNames = get(handles.popUpCoVarModel,'String');
        
%         Squared Exponential
%         Squared Exponential with Nugget
%         p-Exponential Function with separate length scale
%         Squared Exponential Kernel  with separate length scale
%         Matern 3/2
%         Matern 3/2 with separate length scale
%         Matern 5/2
%         Matern 5/2 with separate length scale
        if sum(testIndices==handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.getCovariogramModelChoice)~=1
            str = 'Only the following model are allowed for Matlabs GPR: ';
            for iTestIndice = 1:length(testIndices)
                str = horzcat(str,modelNames{testIndices(iTestIndice)});
                if iTestIndice<length(testIndices)
                    str = horzcat(str,', ');
                end
            end
            errordlg(str)
        else 
            handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.generateRegressionGPModel
        end
    otherwise
end
if exist('mBox1')
    close(mBox1)
end
msgbox('Optimization complete');

[handles] = setCoVarParaTableEntries(handles);

% Delete all callbacks for this button which are made during the
% optimization
drawnow

% Update the GUI and its objects
guidata(hObject, handles);  


% --- Executes on button press in pushbuttonPlotVariogram.
function pushbuttonPlotVariogram_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonPlotVariogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.getDistInput)
    msgbox('(Co-)Variogram estimation has not been excectuted yet')
else
    handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.plotVariogram
    if(handles.KrigingAnalysisObj.KrigingObjects{handles.currentObj}.getnInputVar>2)
        msgbox('Plot Variogram for each Input variable. However, do not interprete to much in these plot since interactions between the input variable are not considered. All distances in the remaining variables are set to zero')
    end
end


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in pushButtonBackButton.
function pushButtonBackButton_Callback(hObject, eventdata, handles)
% hObject    handle to pushButtonBackButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf)


% --- Executes on button press in pushbuttonPlotQuantile.
function pushbuttonPlotQuantile_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonPlotQuantile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.KrigingAnalysisObj.plotQuantilPlot(handles.currentObj)
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
