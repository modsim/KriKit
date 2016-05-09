function [] = calcScreeningAnalysis(obj,varargin)
% [] = calcScreeningAnalysis(obj,KrigingObjectIndex)
%
% Performs screening analysis: Calculation of 3D Kriging estimation for all
% pairwise combinations of input variables. The remaining (nInputvar-2)
% variables are set constant its reference value defined by
% "ReferencePoint".
%
% ReferencePoint ...
% 
% Input:
% - KrigingObjectIndex ... indices of kriging objects of interest [1XnObj]
% 
% Output: -
% 
% You can set: 
% - ReferencePoint ... remaining input variable whihch are not continously
%                      varied over the x- and y-axis are set to values
%                      equal the entries of the "ReferencePoint". (1XnInputVar)
% 
% You can get: -
% 
% For further setting regarding interpolation see documentation of
% "calcMutualInterpolation_23D()".
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.
    
%% Initialization
KrigingObjectIndex = varargin{1};

% Check if everything is correct
checkInput();
    
% calculate pairwise combinations
InputVarCombination = nchoosek(1:obj.KrigingObjects{KrigingObjectIndex}.getnInputVar,2);
nCombination = nchoosek(obj.KrigingObjects{KrigingObjectIndex}.getnInputVar,2);

% Allocate memory
obj.KrigingPrediction_Screening{KrigingObjectIndex,1} = cell(nCombination,1);
obj.KrigingPrediction_Screening{KrigingObjectIndex,2} = cell(nCombination,1);
obj.KrigingPrediction_Screening{KrigingObjectIndex,3} = cell(nCombination,1);
    
%% Do Kriging estimation
for iCombinations = 1:nCombination
    % Define the inputVariable which are not changed continously the on the
    % x and y-axis in current plot
    indicesNonContiInputVar = setdiff(1:obj.KrigingObjects{KrigingObjectIndex}.getnInputVar,InputVarCombination(iCombinations,:));

    % Calculate plot
    obj.calcInterpolation_3D(KrigingObjectIndex,InputVarCombination(iCombinations,:),indicesNonContiInputVar,obj.ReferencePoint(indicesNonContiInputVar));

    % Save result
    obj.KrigingPrediction_Screening{KrigingObjectIndex,1}{iCombinations} = obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,1};
    obj.KrigingPrediction_Screening{KrigingObjectIndex,2}{iCombinations} = obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,2};
    obj.KrigingPrediction_Screening{KrigingObjectIndex,3}{iCombinations} = obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,3};
end

%% Nested Function
function [] = checkInput()

   
    if length(KrigingObjectIndex)>1
        error('Only one Kriging Object can be considered ("KrigingObjectIndex" must of size 1X1)')
    end
    % Determine the input values
    defineBoundOfInputVar(obj,KrigingObjectIndex);
    if length(obj.LBInputVarInterpolation{KrigingObjectIndex})~=obj.KrigingObjects{KrigingObjectIndex}.getnInputVar
        error('"LBInputVarInterpolation has to be a nInputVar x 1 array"')
    end
    if length(obj.UBInputVarInterpolation{KrigingObjectIndex})~=obj.KrigingObjects{KrigingObjectIndex}.getnInputVar
        error('"UBInputVarInterpolation has to be a nInputVar x 1 array"')
    end
    if length(obj.ReferencePoint)~=obj.KrigingObjects{KrigingObjectIndex}.getnInputVar
        error('"ReferencePoint" has be a vector of size nInputVar(=%i)',obj.KrigingObjects{KrigingObjectIndex}.getnInputVar)
    end
    if obj.ShowBasisFct==1
        warning('"ShowBasisFct"=1 and therefore the screening analysis is only done for the basis function')
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
