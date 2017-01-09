function [] = calcInterpolation_nD(obj,varargin)
% calcInterpolation_nD(KrigingObjectIndex,InputVarIndices,ConstantInputVarIndices,ConstantInputVarValue)
%
% Input
% KrigingObjectIndex ... index of the kriging object which should
%                        be used in this function.
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

% InputVarIndices ... Is a Matrix which contains the indices of the input
%  variable which should be analyzed in this interpolation. The rows are
%  associated with the rows in the subplot. The first three columns
%  represent indices of the input variables which are continouesly varied
%  of the interpolation on the x- and y-axis. The fourth column is the
%  index of input variable which is discrete plotted over the columns of
%  the subplot. (nRowsOfSubplotX3)

% ConstantInputVarIndices,ConstantInputVarValue ... if more then 4
%   input variables exist then the remaining input variables have to be set
%   to a defined value. ConstantInputVarIndices is a matrix that contains
%   the indices of theses input variables and ConstantInputVarValue is a
%   matrix that contains its value. The rows are associated with the rows
%   in the subplot. (nRowsX(ninputVar-4))
%
% Note: If you have only 3 input variables there is no need to set one
%   variable fix. Therefore, for nPlots=1 InputVarIndices has only 3
%   columns instead of 4 and ConstantInputVarIndices = [] and
%   ConstantInputVarValue=[]
%
% For further details about the Kriging interpolation see documentation of
% "calcInterpolation()" 
% 
% You can set:- 
% nPlots ... defines the step size between different value of the input
%            variables InputVarIndices(3) and InputVarIndices(4). The plot
%            is a nPlots X nPlots subplot
% 
% 
% You can get: -
% 
% Copyright 2014-2015: Lars Freier, Eric von Lieres
% See the license note at the end of the file.
        
%% Initialization
KrigingObjectIndex = varargin{1};
inputVarIndices = varargin{2};
nRows = size(inputVarIndices,1);
if length(varargin)>=4
    setConstantIndex = varargin{3};
    setConstantValue = varargin{4};
    if ~isempty(setConstantIndex)&&size(setConstantIndex,1)~=nRows&&size(setConstantValue,1)~=nRows
        error('setConstantIndex and setConstantValue should have the same number of rows as inputVarIndices(=%i) but it has %i and %irows',size(setConstantIndex,1),size(setConstantValue,1))
    end
else
    setConstantIndex = [];
    setConstantValue = [];
end
nInputVarConti = 3;

KrigingObjectIndex = obj.checkKrigingIndizes(KrigingObjectIndex);

% Assumption: number of input parameters is constant for all
% chosen input variables
nInputData = obj.KrigingObjects{KrigingObjectIndex(1)}.getnInputVar;
nInputVar = obj.KrigingObjects{KrigingObjectIndex(1)}.getnInputVar;

% Check if enough input variables are defined for the contour plot
if obj.getnPlots>1
    if size(inputVarIndices,2)~=3&&length(varargin)>=3
        error('inputVarIndices must be a matrix of size n(=number of rows in the plot) x 3 but it is  %ix%i ',size(inputVarIndices,1),size(inputVarIndices,2));
    end
end

% Check if all input variables are associated with continous variation or
% staying constant
if (size(inputVarIndices,2)+size(setConstantIndex,2))~=nInputVar
    error('Sum of entries in inputVarIndices and ConstantInputVarIndices has to equal to the number of inputvariable %i but it is (%i+%i)',nInputVar,size(inputVarIndices,2),size(setConstantIndex,2))
end

% Use bounds of input variable for variation over columns
obj.defineBoundOfInputVar(KrigingObjectIndex);
LB = obj.getLBInputVarInterpolation{KrigingObjectIndex};
UB = obj.getUBInputVarInterpolation{KrigingObjectIndex};

% Check if input is ok
for iKriging = KrigingObjectIndex'
    if size(inputVarIndices,2)+size(setConstantIndex,2)~=obj.KrigingObjects{iKriging}.getnInputVar
        error('Entries in setConstantIndex and inputVarIndices must be together equal to the number of input variables but they are %i + %i',size(inputVarIndices,2),size(setConstantIndex,2))
    end
end

if size(setConstantIndex,2)==0&&obj.getnPlots~=1&&size(inputVarIndices,2)~=4
    obj.setnPlots(1);
    warning('No ConstantInputVarIndices and/or ConstantInputVarValue was defined. nPlots is set to 1')
end

InputMatrixProto = cell(nInputData,nRows);
InputMatrix = zeros(obj.getAccuracy^(nInputVarConti-1)*obj.getnPlots*nRows,nInputData);

%% Do Predictions
for iRow = 1:nRows
    
    % Variable with fine grid
    for iInputVar = inputVarIndices(iRow,1:3)
        InputMatrixProto{iInputVar,iRow} = linspace(LB(iInputVar),UB(iInputVar),obj.getAccuracy);
    end

    if length(varargin)>=2
        
        % Variable with coarse grid
        if size(inputVarIndices,2)>3
            for iInputVar = inputVarIndices(iRow,4)
                InputMatrixProto{iInputVar,iRow} = linspace(LB(iInputVar),UB(iInputVar),obj.getnPlots);
            end
        end

        % Variable which stays constant
        for iInputVar = 1:size(setConstantIndex,2)
            InputMatrixProto{setConstantIndex(iRow,iInputVar),iRow} = setConstantValue(iRow,iInputVar);
        end
        
        if ~isempty(find(cellfun(@isempty,InputMatrixProto(:,iRow))==1,1))
            error('You have to define each input variable either as constant(setConstantIndex) or varing(inputVarIndices) for each row in the plot')
        end

    end

    % Actual Interpolation
    for iPlot = 1:obj.getnPlots
        % Positions in the matrix
        startIndex = obj.getAccuracy^(nInputVarConti-1)*obj.getnPlots*(iRow-1) + (iPlot-1)*obj.getAccuracy^(nInputVarConti-1);
        endIndex = obj.getAccuracy^(nInputVarConti-1)*obj.getnPlots*(iRow-1) + iPlot*obj.getAccuracy^(nInputVarConti-1);
        
        % Update Matrix
        InputMatrix(startIndex+1:endIndex,inputVarIndices(4)) = InputMatrixProto{inputVarIndices(4),iRow}(iPlot);
        InputMatrix(startIndex+1:endIndex,setConstantIndex) = setConstantValue;
        
        for iComponent = 1:nInputVarConti-1
            uniqueRows=unique(InputMatrix(startIndex+1:endIndex,:),'rows');
            nUnique = size(uniqueRows,1);
            iIndex = 0;
            for iUnique = 1:nUnique
                concVec = linspace(0,1-sum(uniqueRows(iUnique,obj.PartOfMixture),2),obj.getAccuracy);
                for iAccuary=1:obj.getAccuracy
                    InputMatrix(iIndex+(iAccuary-1)*obj.getAccuracy^(nInputVarConti-1-iComponent)+startIndex + 1:iIndex+(iAccuary)*obj.getAccuracy^(nInputVarConti-1-iComponent)+startIndex,inputVarIndices(iComponent)) = ...
                        concVec(iAccuary);
                end
                iIndex = iIndex+obj.getAccuracy^(iComponent-1);
            end
        end
        
        % Use definition of mixtrue: Sum of input variables has to be 100%
        
        idxAll = 1:obj.KrigingObjects{KrigingObjectIndex}.getnInputVar;
        idxDoNotMod = [inputVarIndices(3),idxAll(~obj.PartOfMixture)];
        remainintVar = setdiff(1:obj.KrigingObjects{KrigingObjectIndex}.getnInputVar,idxDoNotMod);
        InputMatrix(startIndex+1:endIndex,inputVarIndices(3)) = 1- sum(InputMatrix(startIndex+1:endIndex,remainintVar),2);

        threshold = [1-1e-10,1+1e-10];
        if sum(sum(InputMatrix(startIndex+1:endIndex,:),2)<=threshold(1)&sum(InputMatrix(startIndex+1:endIndex,:),2)>=threshold(2))>=1
            warning('Sum of concentration is not 1')
            keyboard
        end
    end
    
end

%% Saving
% Get kriging approximation at these points and save all
% important results
Output = obj.doMututalPrediction(KrigingObjectIndex,InputMatrix);
for iKriging=KrigingObjectIndex'
    obj.KrigingPrediction_InterpolationnD{iKriging,1} = Output;
    obj.KrigingPrediction_InterpolationnD{iKriging,2} = InputMatrix;
    obj.KrigingPrediction_InterpolationnD{iKriging,3} = inputVarIndices;
    obj.KrigingPrediction_InterpolationnD{iKriging,4} = InputMatrixProto;
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
