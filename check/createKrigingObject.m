% clear all
KrigingObject = AnalyzeKriging();
% Copyright 2014-2015: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

KrigingObject.addKrigingObject(2,'TestObj');
KrigingObject.KrigingObjects{1}.setNormInput(1)
KrigingObject.KrigingObjects{1}.setNormOutput(1)
KrigingObject.KrigingObjects{1}.setInputData(load(strcat(currentDir,'/check/testData/input1.txt')))
KrigingObject.KrigingObjects{1}.setOutputData(load(strcat(currentDir,'/check/testData/output1.txt'))')

disp('Kriging Object was created')
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
