function [ ] = solveLeastSquareCovariogramFminSearch(obj)
% Find the optimal parameter Set which for fitting the variogram
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.
    
    % Intitialize Variogram Parameters
    obj.InitializeCovariogramParameters()

    %% Calculate
    options = obj.setOptionsFminSearch();
    obj.covariance_zero = [];
    switch obj.CovariogramEstimationType
        case 1
            [optPara minError] ...
             = fminsearch(@obj.LeastSquareFunctionCovariogram,...
                       obj.InitialCovariogramParameters,...
                       options);
        case 2
            [optPara minError] ...
             = fminsearch(@obj.calcCrossOverQuality,...
                       obj.InitialCovariogramParameters,options);
       case 3
            [optPara minError] ...
             = fminsearch(@obj.MarginalLikelihoodCovarPara,...
                       obj.InitialCovariogramParameters,...
                       options);
        otherwise
            error('CovariogramEstimationType=%d is not defined',obj.CovariogramEstimationType)
    end
    
    % Recalculate Covariagram Matrix if prediction is applied
    obj.checkVariogram = 0;
    % Save Results
    obj.saveOptimizedVariogramParameters(optPara)
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
