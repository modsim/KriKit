function [ ] = solveLeastSquareBasisFct(obj)
% [] = solveLeastSquareBasisFct
% Find the optimal nonlinear basis function parameter by minimizing least
% square error between data and basis function. (Local Optimization via
% fminsearch /fmincon). fminsearch is used when no bounds are defined (see
% below). 
%
% You can set: 
% - BasisFct ... set of basis functions
% - InitialBasisFctParameters ... values for the nonlinear basis fucntion
%       parameter, used local optimizer
% - LBBasisFctParameters/UBBasisFctParameters ... upper and lower bounds of
%                                                 nonlinear basis function
%                                                 parameter
% - Optimizer options see setOptionsFminCon or setOptionsFminSearch
%
% You can get: 
% - ModelErrorBasisFct ... final model error
% - BasisFctParameters ... fina parameter set
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.
    
    %% Calculate
    obj.covariance_zero = [];
    
    % Initialization
    if isempty(obj.InitialBasisFctParameters)
        obj.InitialBasisFctParameters = ones(obj.nBasisFctParameters,1);
    end
    
    % Actual Optimization
    if isempty(obj.LBBasisFctParameters)&&isempty(obj.UBBasisFctParameters)
        options = obj.setOptionsFminSearch();
        	[optPara , minError] ...
                          = fminsearch(@obj.LeastSquareFunctionBasisFct,...
                            obj.InitialBasisFctParameters,...
                            options);
    else
        options = obj.setOptionsFminCon();
        [optPara, minError] ...
         = fmincon(@obj.LeastSquareFunctionBasisFct,...
                   obj.InitialBasisFctParameters,...
                   [],[],[],[],...
                   obj.LBBasisFctParameters,obj.UBBasisFctParameters,...
                   [],options);
    end
    
    obj.setBasisFctParameters(optPara);
    obj.ModelErrorBasisFct = minError;

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
