function  [isInInvestigatedSection]=findIndicesOfDataForPlot(obj,KrigingObjectIndex,Data)
% [isInInvestigatedSection]=findIndicesOfDataForPlot(obj,KrigingObjectIndex,Data)
%
% Determines which data points are in the defined interpolation range
%
% KrigingObjectIndex ... Kriging object which is currently investigated
% Data ... matrix containing input values whcih are to be checked
%          (nTestPointsXnInputVar) 
%
% You can set: -
%
% 
% You can get: -
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

isInInvestigatedSection = true(size(Data,1),1);
for iData = 1:size(Data,1)
    for iVar=1:size(Data,2)
        isInInvestigatedSection(iData) = isInInvestigatedSection(iData)&...
                   Data(iData,iVar)>=obj.LBInputVarInterpolation{KrigingObjectIndex}(iVar) &...
                   Data(iData,iVar)<=obj.UBInputVarInterpolation{KrigingObjectIndex}(iVar);
    end 
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
