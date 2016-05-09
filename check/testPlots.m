nInVar = KrigingObject.KrigingObjects{1}.getnInputVar;
% Copyright 2014-2015: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

KrigingObject.KrigingObjects{1}.setCovariogramModelChoice(4)
for iUseGPR = 1:2
    if iUseGPR==1
        if license('checkout','Statistics_Toolbox');
            KrigingObject.KrigingObjects{1}.setUseMatlabRegressionGP(true)
            KrigingObject.KrigingObjects{1}.generateRegressionGPModel
        else
            continue
        end
    else
        KrigingObject.KrigingObjects{1}.setUseMatlabRegressionGP(false)
    end
    
    KrigingObject.calcAndPlotInterpolation_2D_BestChoice(1,1);
    close all
    KrigingObject.calcAndPlotInterpolation_3D_BestChoice(1,[1,2]);
    close all

    KrigingObject.calcInterpolation_2D(1,1,2:nInVar,ones(1,nInVar-1));
    KrigingObject.plotInterpolation_2D(1)
    close all

    KrigingObject.calcInterpolation_3D(1,[1,2],3:nInVar,ones(1,nInVar-2));
    KrigingObject.plotInterpolation_3D(1)
    close all

    KrigingObject.calcInterpolation_nD(1,[1,2,3],[4,5],[1,1]);
    KrigingObject.plotInterpolation_nD(1)
    close all

    KrigingObject.setReferencePoint(zeros(1,5))
    KrigingObject.calcScreeningAnalysis(1);
    KrigingObject.plotScreeningAnalysisExpectedImprovement(1)
    KrigingObject.plotScreeningAnalysisKrigingInterpolation(1)
    close all
end

disp('Plot test finished')
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
