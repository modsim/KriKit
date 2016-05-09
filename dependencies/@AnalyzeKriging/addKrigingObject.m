function []=addKrigingObject(obj,KrigingType,KrigingObjectName)
% []=AddKrigingObject(KrigingType,KrigingObjectName)
%
% Add new Kriging object. The Object has the type
%
% Input:
% KrigingType ... integer which defines the type of added Kriging Object:
%                 KrigingType=1 ... OrdinaryKriging
%                 KrigingType=2 ... UniversalKriging
%                 KrigingType=3 ... TaylorKriging
% KrigingObjectName ... String Array which defines the name of the
%                       objective (output variable) 
%
% You can set: -
%
% You can get: - 
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.


    switch KrigingType
        case 1
            obj.KrigingObjects{end+1} = OrdinaryKriging;
        case 2
            obj.KrigingObjects{end+1} = UniversalKriging;
        case 3
            obj.KrigingObjects{end+1} = TaylorKriging;
        case 4
            warning('NonstationaryKriging is in the experimental state')
            obj.KrigingObjects{end+1} = NonstationaryKriging;
        otherwise
            error('"KrigingType"=%d was not defined',KrigingType)
    end


    obj.KrigingObjectTypes(end+1) = KrigingType;
    obj.KrigingObjectNames{end+1} = KrigingObjectName;
    obj.InputVarNames{end+1} = {};
    obj.LBInputVarInterpolation{end+1} = [];
    obj.UBInputVarInterpolation{end+1} = [];
    obj.nKrigingObjects = obj.nKrigingObjects + 1;
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
