function [] = defineBoundOfInputVar(obj,varargin)
% [] = defineBoundOfInputVar(obj,KrigingObjectIndex)
% Fill LBInputVarInterpolation and UBInputVarInterpolation is empty
% You can set: -
%
% You can add: - 
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

    KrigingObjectIndex = varargin{1};
    for iKriging=1:length(KrigingObjectIndex)
        % Determin the input values
        if isempty(obj.LBInputVarInterpolation{KrigingObjectIndex(iKriging)})
            obj.LBInputVarInterpolation{KrigingObjectIndex(iKriging)} = min(obj.KrigingObjects{KrigingObjectIndex(iKriging)}.getInputData);
        end
        
        if isempty(obj.UBInputVarInterpolation{KrigingObjectIndex(iKriging)})
            obj.UBInputVarInterpolation{KrigingObjectIndex(iKriging)} = max(obj.KrigingObjects{KrigingObjectIndex(iKriging)}.getInputData);
        end
        
%         if sum(obj.LBInputVarInterpolation{KrigingObjectIndex(iKriging)}>=obj.UBInputVarInterpolation{KrigingObjectIndex(iKriging)})>=1
%             error('One of the upper bounds is smaller or equal to the associated lower bound. ')
%         end
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
