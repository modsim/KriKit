function []=removeKrigingObject(obj,Indices)
% []=RemoveKrigingObject(obj,Indices)
% Removes the KrigingObjects which are identified by its
% Indices in the "KrigingObjects"-array. Be aware that using
% this operation, the indices are decreased of all the Kriging object after 
% the removed object 
%
% Input:
% Indices ... vector of indices of Kriging object which shall be removed.
%             If Indices=[], all KrigingObjects are removed
% 
% You can set: -
%
% You can get: - 
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

    if max(Indices)>length(obj.KrigingObjects)
        error('Maximal Index can be %i',length(obj.KrigingObjects));
    end
    indices=setdiff(1:length(obj.KrigingObjects),Indices);

    % Copy objects
    if ~isempty(indices)
        KrigingObjectsProto = {obj.KrigingObjects{indices}};
        KrigingObjectNamesProto = obj.KrigingObjectNames(indices);
        KrigingObjectTypesProto = obj.KrigingObjectTypes(indices);
        LBInputVarInterpolationProto = {obj.LBInputVarInterpolation{indices}};
        UBInputVarInterpolationProto = {obj.UBInputVarInterpolation{indices}};
        InputVarNamesProto = obj.InputVarNames{indices};

        % Clear and renew objects
        obj.KrigingObjects = KrigingObjectsProto;
        obj.KrigingObjectNames = KrigingObjectNamesProto;
        obj.KrigingObjectTypes = KrigingObjectTypesProto;
        obj.InputVarNames = {InputVarNamesProto};
        obj.LBInputVarInterpolation = LBInputVarInterpolationProto;
        obj.UBInputVarInterpolation = UBInputVarInterpolationProto;

    else
        % Clear and renew objects
        obj.KrigingObjects = {};
        obj.KrigingObjectNames = {};
        obj.KrigingObjectTypes = [];
        obj.LBInputVarInterpolation = cell(0);
        obj.UBInputVarInterpolation = cell(0);
        obj.InputVarNames = cell(0);
    end

    obj.nKrigingObjects = length(obj.KrigingObjects);
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
