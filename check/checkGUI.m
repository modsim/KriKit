% Copyright 2014-2015: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

clc
clear

%% Save Directory
cd ..
currentDir = cd;
cd check\


%% Test if Licenses are available
testLicenses
%
createKrigingObject
%
estimateVariogram
%
testPlots
%
testANOVA
%
clc
disp('Check has finished')

% =============================================================================
%  KriKit - Kriging toolKit
%  
%  Copyright 2014-2015: Lars Freier(1), Eric von Lieres(1)
%                                      
%    (1)Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
%  
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
