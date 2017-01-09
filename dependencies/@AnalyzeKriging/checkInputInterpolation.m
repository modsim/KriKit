function [KrigingObjectIndex,InputVar1,InputVar2,RemainingIndices,RemainingValues] = checkInputInterpolation(obj,varargin)
%[KrigingObjectIndex,InputVar1,InputVar2,RemainingIndices,RemainingValues]=checkInputInterpolation(KrigingObjectIndex,inputVarIndices,RemainingIndices,RemainingValues)
% Input:
% KrigingObjectIndex ... vector containing the indices of the kriging
%                        objects of interest
% inputVarIndices ... vector containing the indices of the input variables
%                     of interest (Continuous variation)
% RemainingIndices ... vector containing the indices of the input variables
%                      hold constant during interpolation
% RemainingValues ... vector containing the values of the input variables
%                      at which they are hold constant during interpolation
% 
% You can set: -
%
% You can get: - 
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.
KrigingObjectIndex = varargin{1};
dimensionInterpolation = varargin{end};
switch obj.KrigingObjects{KrigingObjectIndex(1)}.getnInputVar
    case 1
        if dimensionInterpolation>2
            error('calcInterpolation_3D is only possible if at least 2 input variables exist')
        end
        InputVar1 = 1;
    case 2
        if length(varargin)>2
            InputVar1 = varargin{2}(1);
            InputVar2 = setdiff(1:2,InputVar1);
        else
            InputVar1 = 1;
            InputVar2 = 2;
        end
    otherwise
        switch  dimensionInterpolation
            case 2
                if length(varargin{2})~=1
                    error('Second input must have a length of 1')
                end
            case 3
                if length(varargin{2})~=2
                    error('Second input must have a length of 2')
                end
                InputVar2 = varargin{2}(2);
            otherwise
                error('unexpected dimension of Interpolation')
        end
        
        InputVar1 = varargin{2}(1);
        
end
if dimensionInterpolation==2 
    InputVar2=[];
end

if length(varargin)>3
    RemainingIndices = varargin{3};
    RemainingValues = varargin{4};
    if length(unique([varargin{2},RemainingIndices]))~=obj.KrigingObjects{KrigingObjectIndex(1)}.getnInputVar
        errorCheck();
        error('Entries in "[InputVar1,InputVar2]" and "RemainingIndices" must be unique')
    end
    if dimensionInterpolation>2&&(size(RemainingValues,1)*size(RemainingValues,2))~=obj.KrigingObjects{KrigingObjectIndex(1)}.getnInputVar-2
        errorCheck();
        error('RemainingValues must be a vector with length %i',obj.KrigingObjects{KrigingObjectIndex(1)}.getnInputVar-2)
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
