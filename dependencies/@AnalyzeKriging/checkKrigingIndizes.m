function [indicesOfKrigingObjects] = checkKrigingIndizes(obj,varargin)
%[indicesOfKrigingObjects] = checkKrigingIndizes(indicesOfKrigingObjects)
% 
% You can set: -
%
% You can get: - 
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

indicesOfKrigingObjects = varargin{1};

% Normalize dimensions of indicesOfKrigingObjects
if ~any(size(indicesOfKrigingObjects)==1)
    error('indicesOfKrigingObjects must be an array')
end
if size(indicesOfKrigingObjects,1)==1&&size(indicesOfKrigingObjects,2)>1
    indicesOfKrigingObjects = indicesOfKrigingObjects';
end

% Check Content
indicesOfKrigingObjects = round(indicesOfKrigingObjects);
if sum(indicesOfKrigingObjects<0|indicesOfKrigingObjects>obj.getnKrigingObjects)>length(indicesOfKrigingObjects)
    error('Only indices between 0 and maximal number of Kriging objects (%i) are allowed',obj.getnKrigingObjects)
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
