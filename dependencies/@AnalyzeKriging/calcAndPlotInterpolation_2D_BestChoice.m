function []=calcAndPlotInterpolation_2D_BestChoice(obj,varargin)
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.
 % []=calcAndPlotInterpolation_2D_BestChoice(KrigingObjectIndex,indices)
%
% Kriging interpolation in the 2D space (two input variable + output
% variable). Value of the fixes in input variable are automatically
% adjusted.
% 
% For further information, see documentation of "calcAndPlotInterpolation_23D_BestChoice()" 
%
% Copyright 2014-2015: Lars Freier, Eric von Lieres
% See the license note at the end of the file


calcAndPlotInterpolation_23D_BestChoice(obj,varargin{1:end},2);
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
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
