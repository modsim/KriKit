function [] = estimateBasisFctParametersLocal(obj)
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
                      = fminsearch(@obj.calcKrigingParaEstimationObjFct,...
                        obj.InitialBasisFctParameters,...
                        options);
    else
        options = obj.setOptionsFminCon();
        [optPara, minError] ...
         = fmincon(@obj.calcKrigingParaEstimationObjFct,...
                   obj.InitialBasisFctParameters,...
                   [],[],[],[],...
                   obj.LBBasisFctParameters,obj.UBBasisFctParameters,...
                   [],options);
    end
    obj.BasisFctParameters = optPara;
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
