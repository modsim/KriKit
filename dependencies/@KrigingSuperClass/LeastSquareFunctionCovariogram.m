function [SQ_CoVar] = LeastSquareFunctionCovariogram(obj,varargin)
% [SQ_CoVar] = LeastSquareFunctionCovariogram(parameterSet)
% Calculated the least square error between of the variance estimation and
% variogram function
%
% You can set: -
% - Covariogram model parameters: sigma, sigmaError, p, theta
% - CovariogramModelChoice .. decide for covariogram model
%
% You can get: -
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

    % Define the model parameter values 
    obj.saveOptimizedVariogramParameters(varargin{1});
    
    % Caculate square error
    obj.covariance_zero = ones(size(obj.DistInput,1),1)*obj.CovarModel(zeros(1,size(obj.DistInput,2)),1);
    covariance  = obj.CovarModel(obj.DistInput,1);
    SQ_CoVar    = ((obj.covariance_zero-covariance)-obj.VarianceEstimation).^2;
    SQ_CoVar    = sum(SQ_CoVar);
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
