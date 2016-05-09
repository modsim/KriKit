function labelPlots_23D(obj,varargin)
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

    KrigingObjectIndex = varargin{1};

    
    if isempty(obj.getInputVarNames(KrigingObjectIndex(1)))
        xString = sprintf('Input Variable %i=%g\n',obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,3}(1),max(obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,2}(:,obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,3}(1))) );
        yString = sprintf('Input Variable %i=%g\n',obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,3}(2),max(obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,2}(:,obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,3}(2))) );
        zString = sprintf('Input Variable %i=%g\n',obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,3}(3),max(obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,2}(:,obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,3}(3))) );
    else
        stringInputVar = obj.getInputVarNames(KrigingObjectIndex(1));
        xString = sprintf('%s =%g\n',stringInputVar{obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,3}(1)},max(obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,2}(:,obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,3}(1))) );
        yString = sprintf('%s =%g\n',stringInputVar{obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,3}(2)},max(obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,2}(:,obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,3}(2))) );
        zString = sprintf('%s =%g\n',stringInputVar{obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,3}(3)},max(obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,2}(:,obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,3}(3))) );
    end
    ternlabel(zString,yString,xString);
        
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
