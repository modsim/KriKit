function []=resetKrigingObject(obj,Indices)
% []=resetKrigingObjec(Indices)
% 
% Reset chosen Kriging Objects to their initial state (right after callcing
% addKrigingObject())
%
% Input:
% Indices ... vector of indices of Kriging object which shall be reseted.
%
% You can set: -
%
% You can get: - 
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

    if any(Indices>obj.getnKrigingObjects)
        error('Indices cannot be bigger than numbero f Kriging Objectives (%i)',obj.nKrigingObjects)
    end

    if ~any(size(Indices)==1)
        error('Indices must be a vector')
    end

    % Add New Objects
    nIndices = length(Indices);
    for iIndices=1:nIndices
        addKrigingObject(obj,obj.KrigingObjectTypes(Indices(iIndices)),obj.KrigingObjectNames{Indices(iIndices)})
    end
    
    % Replace old ones 
    obj.KrigingObjects{Indices} = obj.KrigingObjects{end-nIndices+1:end};
    removeKrigingObject(obj,obj.nKrigingObjects-nIndices+1:obj.nKrigingObjects)
    
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
