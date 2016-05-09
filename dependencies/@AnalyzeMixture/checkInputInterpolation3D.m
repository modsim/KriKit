function [KrigingObjectIndex,InputVar1,InputVar2,InputVar3,RemainingIndices,RemainingValues] = checkInputInterpolation3D(obj,varargin)
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.
    KrigingObjectIndex = varargin{1};
    if obj.KrigingObjects{KrigingObjectIndex}.getnInputVar<3
        error('calcInterpolation_3D (mixture) is only possible if at least 3 input variables exist')
    else
        InputVar1 = varargin{2}(1);
        InputVar2 = varargin{2}(2);
        InputVar3 = varargin{2}(3);
    end
    
    if length(varargin)>2
        RemainingIndices = varargin{3};
        RemainingValues = varargin{4};
        if length(unique([varargin{2},RemainingIndices]))~=obj.KrigingObjects{KrigingObjectIndex}.getnInputVar
            error('Entries in "[InputVar1,InputVar2,InputVar3]" and "RemainingIndices" must be unique')
        end
        if (size(RemainingValues,1)*size(RemainingValues,2))~=obj.KrigingObjects{KrigingObjectIndex}.getnInputVar-3
            error('RemainingValues must be a vector with length %i',obj.KrigingObjects{KrigingObjectIndex}.getnInputVar-3)
        end

        % Transpose if neccesary
        if size(RemainingIndices,1)>size(RemainingIndices,2)
            RemainingIndices = RemainingIndices';
        end

        if size(RemainingIndices,1)>size(RemainingIndices,2)
            RemainingIndices = RemainingIndices';
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
