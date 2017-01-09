function []=setLBUBInputVarInterpolation(obj,KrigingObjectIndex,LBUBInputVarInterpolation,LBorUB)
% [inputData,indexValid] = doTestForOptimality(obj,KrigingIndex,dimensionInterpolation,testValue)
% 
% 
% Input:
%
% You can set: -
%
% You can get: -
% 
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.
    
    % Check Input
    if length(KrigingObjectIndex)~=1||KrigingObjectIndex<1||KrigingObjectIndex>obj.getnKrigingObjects
        error('KrigingObjectIndex must be a scalar between 1 and %i',obj.getnKrigingObjects)
    end
    
    % A row vector is needed
    if size(LBUBInputVarInterpolation,1)==obj.KrigingObjects{KrigingObjectIndex}.getnInputVar&&size(LBUBInputVarInterpolation,2)==1
        LBUBInputVarInterpolation = LBUBInputVarInterpolation';
    elseif size(LBUBInputVarInterpolation,2)==obj.KrigingObjects{KrigingObjectIndex}.getnInputVar&&size(LBUBInputVarInterpolation,1)==1

    else
        error('LB/UBInputVarInterpolation has to be of dimension 1x%i',obj.KrigingObjects{KrigingObjectIndex}.getnInputVar)
    end

    

    switch LBorUB
        case -1
            obj.LBInputVarInterpolation{KrigingObjectIndex} = LBUBInputVarInterpolation;
        case 1
            if any(obj.LBInputVarInterpolation{KrigingObjectIndex}>LBUBInputVarInterpolation)
                error('One of the upper bounds is smaller than the associated lower bound. Please define the lower bound first')
            end
            if all(obj.LBInputVarInterpolation{KrigingObjectIndex}==LBUBInputVarInterpolation)
                error('All upper bounds is smaller or equal to the associated lower bound. Please define the lower bound first')
            end
            obj.UBInputVarInterpolation{KrigingObjectIndex} = LBUBInputVarInterpolation;
        otherwise
            error('LBorUB has to be -1 (LB) or +1 (UB)')
    end
    

% % %     % This is need for all data set which use an old version of
% % %     % KriKit
% % %     if strcmp(ex.identifier,'MATLAB:cellAssToNonCell')
% % %         obj.LBUBInputVarInterpolation = cell(obj.getnKrigingObjects,1);
% % %         obj.LBUBInputVarInterpolation{KrigingObjectIndex} = LBUBInputVarInterpolation;
% % %     else
% % %         error('Unkown Error during setting LBInputVarInterpolation')
% % %     end


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
