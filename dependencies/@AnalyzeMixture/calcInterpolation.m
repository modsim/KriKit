function []=calcInterpolation(obj,varargin)
% [] = calcInterpolation(KrigingObjectIndex,inputVariableIndices,RemainingIndices,RemainingValues)
% 
% Kriging interpolation of the given gaussian process defined by the
% Kriging objects with the indices KrigingObjectIndex. Interpolation is for
% varying continuously the values of three input variables. Remaining input variables
% are set constant to defined values.
%
% Input: 
% KrigingObjectIndex ... index of the Kriging object which are used. If
%                        "KrigingObjectIndex" is a vector, mutual Kriging
%                        is performed. For more details see documentation
%                        of "doMututalPrediction()"
% inputVariableIndices ... indices of the input variables which are
%                          continuously varied (1X3)
% RemainingIndices ... Array of the remaining input variables. (1X(nInputvar-3))
% RemainingValues ... vector which contains constant input values of
%                     remaining input variables (defined by
%                     RemainingIndices). (1X(nInputvar-3))
% 
% NOTE: Order of RemainingIndices and RemainingValues must be consistent.
% Each entry in RemainingValues is associated with the input variable with
% the index defined by RemainingIndices at the same position
%
% You can set:
% - Accuracy ... defines fineness of grid of input variables used for
%                interpolation. Represents number of steps for linearly
%                spacing between LBInputVarInterpolation and
%                UBInputVarInterpolation
% - ShowBasisFct ... if true, then the trend function is plotted
%                    instead of the KrigingPrediction
% 
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

%% Initialization
[KrigingObjectIndex,InputVar1,InputVar2,InputVar3,RemainingIndices,RemainingValues] = checkInputInterpolation3D(obj,varargin{1:end});
inputVarIndices = [InputVar1,InputVar2,InputVar3];
defineBoundOfInputVar(obj,KrigingObjectIndex);

Input = zeros(obj.getAccuracy^2,obj.KrigingObjects{KrigingObjectIndex}.getnInputVar);
Input(:,RemainingIndices) = RemainingValues;

%% Generate Input Grid 
for iComponent = 1:2
    uniqueRows=unique(Input,'rows');
    nUnique = size(uniqueRows,1);
    iIndex = 0;
    for iUnique = 1:nUnique
        concVec = linspace(0,1-sum(uniqueRows(iUnique,obj.PartOfMixture),2),obj.getAccuracy);
        for iAccuary=1:obj.getAccuracy
            indexStart = iIndex+(iAccuary-1)*obj.getAccuracy^(2-iComponent) + 1;
            indexEnd = iIndex+(iAccuary)*obj.getAccuracy^(2-iComponent);
            Input(indexStart:indexEnd,inputVarIndices(iComponent)) = concVec(iAccuary);
        end
        iIndex = iIndex+obj.getAccuracy^(iComponent-1);
    end
end
% Input = createNDGRID([0,0],[1,1],obj.getAccuracy);

% Use definition of mixtrue: Sum of input variables has to be 100%
idxAll = 1:obj.KrigingObjects{KrigingObjectIndex}.getnInputVar;
idxDoNotMod = [inputVarIndices(3),idxAll(~obj.PartOfMixture)];
indicesOthers = setdiff(1:obj.KrigingObjects{KrigingObjectIndex}.getnInputVar,idxDoNotMod);
Input(:,inputVarIndices(3)) = 1- sum(Input(:,indicesOthers),2);
Input = unique(Input,'rows');

%% Make Prediction and save Results
if obj.getShowBasisFct==1
%     BasisFct = obj.KrigingObjects{KrigingObjectIndex}.getBasisFct;
%     if length(BasisFct)>1
%         error('More than one basis function is defined')
%     end
%     
%     % Every basis function is weighted with a Kriging calculated
%     % coefficient
%     obj.KrigingObjects{KrigingObjectIndex}.estimateBasisFctCoefficients;
%     obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,1} = obj.KrigingObjects{KrigingObjectIndex}.getBasisFctCoefficients*...
%                                                                   [BasisFct{1}(obj.KrigingObjects{KrigingObjectIndex}.getBasisFctParameters,Input),...
%                                                                    zeros(size(Input,1),1)];
                                                               
    evaluateBasisFunction(obj,KrigingObjectIndex,Input,size(Input,2));
    
else
    
obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,1} = doMututalPrediction(obj,KrigingObjectIndex,Input);
    
end
obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,2} = Input;
obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,3} = [InputVar1,InputVar2,InputVar3];

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
