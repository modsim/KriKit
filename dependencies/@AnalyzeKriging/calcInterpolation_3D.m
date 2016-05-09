function []=calcInterpolation_3D(obj,varargin)
% [] = calcInterpolation_3D(KrigingObjectIndex,[InputVar1,InputVar2],RemainingIndices,RemainingValues)
%
% For Further Detail see documentation of "calcMutualInterpolation_23D()"
% You can set: -
%
% You can get: - 
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

dimension = 3;
[KrigingObjectIndex,InputVar1,InputVar2,RemainingIndices,RemainingValues] = checkInputInterpolation(obj,varargin{1:end},dimension);
calcMutualInterpolation_23D(obj,KrigingObjectIndex,[InputVar1,InputVar2],RemainingIndices,RemainingValues,dimension);
    
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
