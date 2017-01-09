function []=calcMutualInterpolation_23D(obj,varargin)
% [] = calcMutualInterpolation_23D(KrigingObjectIndex,inputVariableIndices,RemainingIndices,RemainingValues,dimensionInterpolation)
%
% Kriging interpolation of the given gaussian process defined by the
% Kriging objects with the indices KrigingObjectIndex. Interpolation is
% for varying continuously the values of one or two input variables, in the
% case of 2D or 3D interpolation, respectively. Remaining input variables
% are set constant to defined values.
%
% Input:
% KrigingObjectIndex ... index of the Kriging object which are used. If
%                        "KrigingObjectIndex" is a vector, mutual Kriging
%                        is performed. For more details see documentation
%                        of "doMututalPrediction()"
% inputVariableIndices ... indices of the input variables which are
%                          continuously varied (1X2, in case of 3D plot
%                          otherwise 1X1) 
% RemainingIndices ... Array of the remaining input variables. (1X(nInputvar-2(1)))
% RemainingValues ... vector which contains constant input values of
%                     remaining input variables (defined by
%                     RemainingIndices). (1X(nInputvar-2(1)))
% dimensionInterpolation ... decide for 2D or 3D case (dimensionInterpolation={2,3})
%
% NOTE: Order of RemainingIndices and RemainingValues must be consistent.
% Each entry in RemainingValues is associated with the input variable with
% the index defined by RemainingIndices at the same position
%
% You can set:
% - LBInputVarInterpolation, UBInputVarInterpolation ... Defines the
%       range of the input variables in which interpolation is performed
% - Accuracy ... defines fineness of grid of input variables used for
%                interpolation. Represents number of steps for linearly
%                spacing between LBInputVarInterpolation and
%                UBInputVarInterpolation
% - ShowBasisFct ... if true, then the trend function is plotted
%                    instead of the KrigingPrediction
%
% You can get: 
% - KrigingPrediction_Interpolation2D/3D ... save the interpolation results
%                                            for more information see
%                                            documentation of
%                                            getKrigingPrediction 
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.
    
%% Initialization
dimensionInterpolation = varargin{5};
KrigingObjectIndex = varargin{1};

[KrigingObjectIndex,InputVar1,InputVar2,RemainingIndices,RemainingValues]=checkInput(KrigingObjectIndex);
nKrigingObjectIndex = length(KrigingObjectIndex);

% Determine the input values
defineBoundOfInputVar(obj,KrigingObjectIndex,RemainingIndices,RemainingValues);

% Assumption: each Kriging objects has the same number of
% input variables. If not it is later detected
nInputData = obj.KrigingObjects{KrigingObjectIndex(1)}.getnInputVar;
LB_Backup=obj.getLBInputVarInterpolation;
UB_Backup=obj.getUBInputVarInterpolation;

[Input] = defineGridofInputVariables(RemainingIndices,RemainingValues);

%% Actual Prediction
doInterpolation();

%% Save Important information
saveInterpolationinformation();

%% Nested Functions
% -------------------------------------------------------------------------
function [KrigingObjectIndex,InputVar1,InputVar2,RemainingIndices,RemainingValues]=checkInput(KrigingObjectIndex)
    
    % Prelimiary Check
    KrigingObjectIndex = obj.checkKrigingIndizes(KrigingObjectIndex);
    
    % Assumption: each Kriging objects has the same number of
    % input variables. If not it is later detected
    nInputData = obj.KrigingObjects{KrigingObjectIndex(1)}.getnInputVar;

    if nInputData==1
        InputVar1 = 1;
        if dimensionInterpolation==3
            error('calcInterpolation_3D is only possible if at least 2 input variables exist')
        end
    elseif nInputData==1&&dimensionInterpolation<3
        InputVar1 = 1;
        InputVar2 = 2;
    else
        if length(varargin{2})~=dimensionInterpolation-1
            error('Second input must have a length of %i',dimensionInterpolation-1)
        end
        InputVar1 = varargin{2}(1);
        if dimensionInterpolation==3
            InputVar2 = varargin{2}(2);
        end
    end

    if dimensionInterpolation==2
        InputVar2 = [];
    end
    
    
    if length(varargin)>2&&~isempty(varargin{3})&&~isempty(varargin{4})
        RemainingIndices = varargin{3};
        RemainingValues = varargin{4};
        if length(unique([varargin{2},RemainingIndices]))~=nInputData
            errorCheck();
            error('Entries in "[InputVar1,InputVar2]" and "RemainingIndices" must be unique')
        end
        if (size(RemainingValues,1)*size(RemainingValues,2))~=nInputData-(dimensionInterpolation-1)
            errorCheck();
            error('RemainingValues must be a vector with length %i',nInputData-(dimensionInterpolation-1))
        end

        % Transpose if neccesary
        if size(RemainingIndices,1)>size(RemainingIndices,2)
            RemainingIndices = RemainingIndices';
        end

        if size(RemainingIndices,1)>size(RemainingIndices,2)
            RemainingIndices = RemainingIndices';
        end
    else
        RemainingIndices = [];
        RemainingValues = [];
    end
%     % Determin the input values
%     defineBoundOfInputVar(obj,KrigingObjectIndex);

end
% -------------------------------------------------------------------------
function [Input] = defineGridofInputVariables(RemainingIndices,RemainingValues)
    
    % First Linear Spacing between lower and upper bound
    if nKrigingObjectIndex==1
        I1=linspace(obj.LBInputVarInterpolation{KrigingObjectIndex}(InputVar1),obj.UBInputVarInterpolation{KrigingObjectIndex}(InputVar1),obj.Accuracy);
        if dimensionInterpolation==3
            I2=linspace(obj.LBInputVarInterpolation{KrigingObjectIndex}(InputVar2),obj.UBInputVarInterpolation{KrigingObjectIndex}(InputVar2),obj.Accuracy);
%             [a,b]=ndgrid(I1,I2);
        end
    else
        LB = [obj.LBInputVarInterpolation{KrigingObjectIndex}];
        if length(LB)~=nKrigingObjectIndex*nInputData
            error('Number of input variables must be constant over all chosen Kriging objects')
        end
        LB = reshape(LB,nInputData,nKrigingObjectIndex)';
        UB = reshape([obj.UBInputVarInterpolation{KrigingObjectIndex}],nInputData,nKrigingObjectIndex)';
        I1=linspaceNDim(min(LB(:,InputVar1)),max(UB(:,InputVar1)),obj.Accuracy);
        if dimensionInterpolation==3
            I2=linspaceNDim(min(LB(:,InputVar2)),max(UB(:,InputVar2)),obj.Accuracy);
%             [a,b]=ndgrid(I1,I2);
        end
    end
    
    % Determine Grid over both input variables
    if dimensionInterpolation==3
        [a,b]=ndgrid(I1,I2);
    end
    

    % Collect all input Data (used in the next step)
    nData = 0;
    for iKrigingNested=KrigingObjectIndex'
        nData = nData+obj.KrigingObjects{iKrigingNested}.getnExperiments;
    end
    InputData = ones(nData,nInputData);
    index=1;
    for iKrigingNested=KrigingObjectIndex'
        InputData(index:index+obj.KrigingObjects{iKrigingNested}.getnExperiments-1,:) = obj.KrigingObjects{iKrigingNested}.getInputData;
    end
    
    % Initialially, define each value of the remaining variable by its
    % median oer all sample points and replace it with user defined value,
    % if given
    Input = bsxfun(@times,ones(obj.Accuracy^(dimensionInterpolation-1),size(InputData,2)),median(InputData));
    switch dimensionInterpolation
        case 2
            Input(:,InputVar1)=I1;
        case 3
            Input(:,InputVar1)=a(:);
            Input(:,InputVar2)=b(:);
        otherwise
            error('Diminsion is unusual')
    end

    % Replace with user defined value
    if length(varargin)>2&&~isempty(varargin{3})&&~isempty(varargin{4})
        Input(:,RemainingIndices)=bsxfun(@times,ones(size(Input(:,RemainingIndices))),RemainingValues);
    end
end
% -------------------------------------------------------------------------
function [] = doInterpolation()
    if obj.ShowBasisFct==1
        % Calculated only evaluation of the basis function
        evaluateBasisFunction(obj,KrigingObjectIndex,Input,dimensionInterpolation);
    else
        % Actual Kriging Estimation
        if nKrigingObjectIndex==1
            OutputMatrix = obj.KrigingObjects{KrigingObjectIndex}.prediction(Input);
        else
            OutputMatrix = obj.doMututalPrediction(KrigingObjectIndex,Input);
        end

        for iKrigingNested=KrigingObjectIndex'
            switch dimensionInterpolation
                case 2
                    obj.KrigingPrediction_Interpolation2D{iKrigingNested,1} = OutputMatrix;
                case 3
                    obj.KrigingPrediction_Interpolation3D{iKrigingNested,1} = OutputMatrix;
                otherwise
                    error('Diminsion is unusual')
            end
        end
    end
end
% -------------------------------------------------------------------------
function [] = saveInterpolationinformation()
    for iKrigingNested=KrigingObjectIndex'

        switch dimensionInterpolation
            case 2
                obj.KrigingPrediction_Interpolation2D{iKrigingNested,2} = Input;
                obj.KrigingPrediction_Interpolation2D{iKrigingNested,3} = InputVar1;
            case 3
                obj.KrigingPrediction_Interpolation3D{iKrigingNested,2} = Input;
                obj.KrigingPrediction_Interpolation3D{iKrigingNested,3} = [InputVar1,InputVar2];
            otherwise
                error('Diminsion is unusual')
        end

        if ~isempty(LB_Backup{iKrigingNested})
            obj.setLBInputVarInterpolation(iKrigingNested,LB_Backup{iKrigingNested});
            obj.setUBInputVarInterpolation(iKrigingNested,UB_Backup{iKrigingNested});
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
