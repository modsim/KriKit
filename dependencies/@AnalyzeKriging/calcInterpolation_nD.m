function [] = calcInterpolation_nD(obj,varargin)
% calcInterpolation_nD(KrigingObjectIndex,InputVarIndices,ConstantInputVarIndices,ConstantInputVarValue)
%
% Input
% KrigingObjectIndex ... index of the kriging object which should
%                        be used in this function.
% InputVarIndices ... Is a Matrix which contains the indices of the input variable
%   which should be analyzed in this interpolation. The rows are associated
%   with the rows in the subplot. The first two columns represent indices
%   of the input variables which are continouesly varied of the
%   interpolation on the x- and y-axis. The third column is the index of
%   input variable which is discrete plotted over the columns of the
%   subplot. (nRowsOfSubplotX3)
% ConstantInputVarIndices,ConstantInputVarValue ... if more then 3
%   input variables exist then the remaining input variables have to be set
%   to a defined value. ConstantInputVarIndices is a matrix that contains
%   the indices of theses input variables and ConstantInputVarValue is a
%   matrix that contains its value. The rows are associated with the rows
%   in the subplot. (nRowsX(ninputVar-3))
%
% Note: If you have only 2 input variables there is no need to set one
%   variable fix. Therefore, for nPlots=1 InputVarIndices has only 2
%   columns instead of 3 and ConstantInputVarIndices = [] and
%   ConstantInputVarValue=[]
%
% For further details about the Kriging interpolation see documentation of
% "calcMutualInterpolation_23D()" 
% 
% You can set:- 
% nPlots ... defines the step size between different value of the input
%            variables InputVarIndices(3) and InputVarIndices(4). The plot
%            is a nPlots X nPlots subplot
% 
% 
% You can get: -
% 
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.


        
%% Initialization
[KrigingObjectIndex,...
 inputVarIndices,...
 nRows,...
 setConstantIndex,...
 setConstantValue,...
 LB,...
 UB,...
 InputMatrixProto,...
 InputMatrix] = doInitialization();

%%
[InputMatrixProto,InputMatrix] = calculateInputMatrix(InputMatrixProto,InputMatrix);

%% Actual Kriging Interpolation
Output = obj.doMututalPrediction(KrigingObjectIndex,InputMatrix);

%% Save Important Information
for iKriging=KrigingObjectIndex'
    obj.KrigingPrediction_InterpolationnD{iKriging,1} = Output;
    obj.KrigingPrediction_InterpolationnD{iKriging,2} = InputMatrix;
    obj.KrigingPrediction_InterpolationnD{iKriging,3} = inputVarIndices;
    obj.KrigingPrediction_InterpolationnD{iKriging,4} = InputMatrixProto;
end

%% Nested Function
% -------------------------------------------------------------------------
function [KrigingObjectIndex,inputVarIndices,nRows,setConstantIndex,setConstantValue,LB,UB,InputMatrixProto,InputMatrix]=doInitialization()
    
    % Check Input
    KrigingObjectIndex = varargin{1};
    KrigingObjectIndex = obj.checkKrigingIndizes(KrigingObjectIndex);
    inputVarIndices = varargin{2};
    nRows = size(inputVarIndices,1);

    % Assumption: number of input parameters is constant for all
    % chosen input variables
    nInputVar = obj.KrigingObjects{KrigingObjectIndex(1)}.getnInputVar;
    
    % Check if enough information are provided
    if obj.nPlots>1
        if size(inputVarIndices,2)~=3&&length(varargin)>=3
            error('inputVarIndices must be a matrix of size n(=number of rows in the plot) x 3 but it is  %ix%i ',size(inputVarIndices,1),size(inputVarIndices,2));
        end
    else
        if size(inputVarIndices,2)~=2
            error('For nPlots==1, inputVarIndices must be a matrix of size n(=number of rows in the plot) x 2 but it is  %ix%i ',size(inputVarIndices,1),size(inputVarIndices,2));
        end
    end
    
    % Define which input variable are kept constanst over the entire row
    % and at which value
    if length(varargin)>=3
        setConstantIndex = varargin{3};
        setConstantValue = varargin{4};
        if ~isempty(setConstantIndex)&&size(setConstantIndex,1)~=nRows&&size(setConstantValue,1)~=nRows
            error('setConstantIndex and setConstantValue should have the same number of rows as inputVarIndices(=%i) but it has %i and %irows',size(setConstantIndex,1),size(setConstantValue,1))
        end
    else
        setConstantIndex = [];
        setConstantValue = [];
    end
        % Check Results
    if (size(inputVarIndices,2)+size(setConstantIndex,2))~=nInputVar
        error('Sum of entries in inputVarIndices and ConstantInputVarIndices has to be equal to the number of inputvariable %i but it is (%i+%i)',nInputVar,size(inputVarIndices,2),size(setConstantIndex,2))
    end

    % Define Lower and upper bound of interpolation range
    obj.defineBoundOfInputVar(KrigingObjectIndex);
    LB = obj.LBInputVarInterpolation{KrigingObjectIndex};
    UB = obj.UBInputVarInterpolation{KrigingObjectIndex};

    % Check if input is ok for all Kriging objects
    for iKrigingNested = KrigingObjectIndex'
        if size(inputVarIndices,2)+size(setConstantIndex,2)~=obj.KrigingObjects{iKrigingNested}.getnInputVar
            error('Entries in setConstantIndex and inputVarIndices must be together equal to the number of input variables but they are %i + %i',size(inputVarIndices,2),size(setConstantIndex,2))
        end
    end

    if size(setConstantIndex,2)==0&&obj.nPlots~=1&&size(inputVarIndices,2)~=3
        obj.nPlots=1;
        warning('No ConstantInputVarIndices and/or ConstantInputVarValue was defined. nPlots is set to 1')
    end
    
    % Allocated Memory
    InputMatrixProto = cell(nInputVar,nRows);
    InputMatrix = zeros(obj.Accuracy^2*obj.nPlots*nRows,nInputVar);
end
% -------------------------------------------------------------------------
function [InputMatrixProto,InputMatrix] = calculateInputMatrix(InputMatrixProto,InputMatrix)
    for iRow = 1:nRows
        % Variable with fine grid
        for iInputVar = inputVarIndices(iRow,1:2)
            InputMatrixProto{iInputVar,iRow} = linspace(LB(iInputVar),UB(iInputVar),obj.Accuracy);
        end
        
        % If variable is defined which is discrete varied of the columns
        nInputVar = obj.KrigingObjects{KrigingObjectIndex(1)}.getnInputVar;
        if length(varargin)>=3||nInputVar==3
            
            % Variable with coarse grid (discrete variation)
            if size(inputVarIndices,2)>2
                for iInputVar = inputVarIndices(iRow,3)
                    InputMatrixProto{iInputVar,iRow} = linspace(LB(iInputVar),UB(iInputVar),obj.nPlots);
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
        
        % Combine Results in InputMatrix
        for iPlot = 1: obj.nPlots
            indexStart = obj.Accuracy^2*obj.nPlots*(iRow-1)+obj.Accuracy^2*(iPlot-1)+1;
            indexEnd = obj.Accuracy^2*obj.nPlots*(iRow-1)+obj.Accuracy^2*iPlot;
            InputMatrix(indexStart:indexEnd,inputVarIndices(iRow,1:2)) = createNDGRID(LB(inputVarIndices(iRow,1:2)),UB(inputVarIndices(iRow,1:2)),obj.Accuracy);
            if obj.nPlots>1
                InputMatrix(indexStart:indexEnd,inputVarIndices(iRow,3)) = InputMatrixProto{inputVarIndices(iRow,3),iRow}(iPlot);
            end
            if ~isempty(setConstantValue)
                InputMatrix(indexStart:indexEnd,setConstantIndex(iRow,:)) = repmat(setConstantValue(iRow,:),obj.Accuracy^2,1);
            end
        end
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
