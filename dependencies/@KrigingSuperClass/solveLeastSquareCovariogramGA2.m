function [ ] = solveLeastSquareCovariogramGA2(obj)
% [] = solveLeastSquareCovariogramGA2
%
% Find the optimal nonlinear basis function parameter by minimizing least
% square error between data and basis function. (Global Optimization via
% Open Source: Genetic Algorithm. Can be downloaded under
% http://www.shef.ac.uk/acse/research/ecrg/gat.html and
% http://codem.group.shef.ac.uk/index.php/ga-toolbox Introcution Paper:
% Chipperfield, A. J., and P. J. Fleming. "The MATLAB genetic algorithm
% toolbox." Applied Control Techniques Using MATLAB, IEE Colloquium on.
% IET, 1995.
%
% You can set: 
% - CovariogramModelChoice ... decide which covariogram model shold be
%                              used. For further information see
%                              documentation of "CovariogramModelChoice" in
%                              "KrigignSuperClass"
% - LBCovariogramModelParameters/UBCovariogramModelParameters... 
%           upper and lower bounds of nonlinear basis function parameter.
%           
% - CovariogramEstimationType ... 
%   1 minimize least squared error between variogram model evaluation and
%     data based variance estimation. See
%     "LeastSquareFunctionCovariogram()" 
%   2 cross validation error. see "calcCrossOverQuality()"
%   3 Maximum likelihood apporach, see "MarginalLikelihoodCovarPara()".
%     This is default.
% - GA-Optimization parameters, see documentation of setOptionsGA
% 
% You can get: 
% - CovariogramModelParameters ... estimated covariogram model parameters
% - ModelErrorVar  ... final model error after optimization
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

    % Intitialize Variogram Parameters
    obj.InitializeCovariogramParameters()

    % Adjust the options
    options=obj.setOptionsGA();
    switch obj.CovariogramEstimationType
        case 1
            obj.covariance_zero = [];      
            [optPara minError] ...
             = obj.ga2(@obj.LeastSquareFunctionCovariogram,length(obj.LBCovariogramModelParameters),...
                       obj.LBCovariogramModelParameters, obj.UBCovariogramModelParameters,...
                       options);
        case 2
            obj.covariance_zero = [];      
            [optPara minError] ...
             = obj.ga2(@obj.calcCrossOverQuality,length(obj.LBCovariogramModelParameters),...
                       obj.LBCovariogramModelParameters, obj.UBCovariogramModelParameters,...
                       options);
        case 3
            obj.covariance_zero = [];    
            [optPara minError] ...
             = obj.ga2(@obj.MarginalLikelihoodCovarPara,length(obj.LBCovariogramModelParameters),...
                       obj.LBCovariogramModelParameters, obj.UBCovariogramModelParameters,...
                       options);
    end
    
    % Recalculate Covariagram Matrix if prediction is applied
    obj.checkVariogram = 0;
    % Save Results
    obj.saveOptimizedVariogramParameters(optPara(1,:))
    obj.ModelErrorVar = minError;

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
