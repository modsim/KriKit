function [inputData,indexValid]=plotOptimum2D(obj,varargin)
% [inputData,indexValid] = plotPlateau2D(obj,KrigingIndex,testValue)
% 
% For further details, see documentation of "plotOptimum23D()"
% 
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

if length(varargin)>1
    testValue = varargin{2};
else
    testValue = [];
end
[inputData,indexValid] = obj.plotOptimum23D(varargin{1},2,testValue);

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
