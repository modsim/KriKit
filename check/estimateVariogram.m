KrigingObject.setPopulationSize(1e1)
KrigingObject.setGenerations(3)
KrigingObject.setTimeLimit(1e3)
KrigingObject.setShowDetails(true)
KrigingObject.KrigingObjects{1}.setCovariogramModelChoice(4)
KrigingObject.KrigingObjects{1}.setBasisFct('polynomial',0)
% Copyright 2014-2015: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

nSolver = 5;
nCovariogramChoices = 8;
nVarEstiType = 3;
solverVec = sort(1:nSolver,'descend');
% warning('off')
for iVar = 1:nVarEstiType
    for iSolver = solverVec
        for iCov = 1:nCovariogramChoices
            
            % Check if license for Optimization Toolbox is available
            if (iSolver==2&&~(license('checkout','Optimization_Toolbox')))||(iSolver==4&&~(license('checkout','GADS_Toolbox')))
                continue
            end
            
            % Set Properties
            KrigingObject.KrigingObjects{1}.setCovariogramEstimationType(iVar)
            KrigingObject.KrigingObjects{1}.setCovariogramModelChoice(iCov)
            KrigingObject.setUseSolver(iSolver);
            
            % Intial Parameters are neccessary for local optimizer
            KrigingObject.KrigingObjects{1}.setInitialCovariogramParameters(KrigingObject.KrigingObjects{1}.getCovariogramModelParameters);
            
            % Set UB sicne otherwise numerical problem may occur
            nPara = length(KrigingObject.KrigingObjects{1}.getUBCovariogramModelParameters);
            KrigingObject.KrigingObjects{1}.setUBCovariogramModelParameters(ones(1,nPara)*2);
            
            
            KrigingObject.KrigingObjects{1}.estimateVariance();
            KrigingObject.KrigingObjects{1}.setnIterationsSolver(10);
            switch iSolver
                case 1
                    KrigingObject.KrigingObjects{1}.solveLeastSquareCovariogramFminSearch;
                case 2
                    if license('checkout','Optimization_Toolbox')
                        KrigingObject.KrigingObjects{1}.solveLeastSquareCovariogram;
                    end
                case 3
                    KrigingObject.KrigingObjects{1}.solveLeastSquareCovariogramGA2;
                case 4
                    if license('checkout','GADS_Toolbox')
                        KrigingObject.KrigingObjects{1}.solveLeastSquareCovariogram;
                    end
                case 5
                    if sum([2,4:8]==iCov)==1&&license('checkout','Statistics_Toolbox');
                        KrigingObject.KrigingObjects{1}.generateRegressionGPModel;
                    end
            end
            
            
        end
    end
end


disp('Covariogram estimation ready')
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
