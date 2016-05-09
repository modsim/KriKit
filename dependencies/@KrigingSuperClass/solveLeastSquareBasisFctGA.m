function [ ] = solveLeastSquareBasisFctGA(obj)
% [] = solveLeastSquareBasisFctGA
% Find the optimal nonlinear basis function parameter by minimizing least
% square error between data and basis function. (Global Optimization via
% mathworks genetic algorithm) 
%
% You can set: 
% - BasisFct ... set of basis functions
% - InitialBasisFctParameters ... values for the nonlinear basis fucntion
%       parameter, used local optimizer
% - LBBasisFctParameters/UBBasisFctParameters ... upper and lower bounds of
%                                                 nonlinear basis function
%                                                 parameter
% - GA-Optimization parameters, see documentation of setOptionsGA
%
% You can get: 
% - ModelErrorBasisFct ... final model error
% - BasisFctParameters ... fina parameter set
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.
    
    % Adjust the options
    options=obj.setOptionsGA();
    
    % Initialization
    if isempty(obj.InitialBasisFctParameters)
        obj.InitialBasisFctParameters = ones(obj.nBasisFctParameters,1);
    end

    obj.covariance_zero = [];   
    [optPara minError] ...
         = ga(@obj.LeastSquareFunctionBasisFct,obj.nBasisFctParameters,...
                   [],[],[],[],...
                   obj.LBBasisFctParameters,obj.UBBasisFctParameters,...
                   [],options);
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
