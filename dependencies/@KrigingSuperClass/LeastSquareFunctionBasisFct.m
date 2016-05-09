function [SQ_CoVar] = LeastSquareFunctionBasisFct(obj,varargin)
% [SQ_CoVar] = LeastSquareFunctionBasisFct(parameterSet)
% Calculated the least square error of between output data and
% basis function estimation
% 
% Note: Linear Parameter are always automatically estimated in the Kriging
% framework. Only Define nonlinear parameter when you "setBasisFct()"
%
% E.g.: use "@(p,x)x(:,1)^p(1)" instead of "@(p,x)p(2)*x(:,1)^p(1)" as the
% function is linear w.r.t. p(2)
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

    % Define the model parameter values 
    switch obj.BasisFctType
        case obj.allowed{1}
            warning('For %s the function LeastSquareFunctionBasisFct will do nothing (no nonlinearparameter)',obj.allowed{1});
        case obj.allowed{end}
            
            if obj.nBasisFctParameters>0
                % Estimate linear parameter via Kriging. Turn off warning
                % for matrix invertations
                warningIDBefore = warning('query','MATLAB:nearlySingularMatrix');
                warning('off','MATLAB:nearlySingularMatrix');
                try
                    obj.estimateBasisFctCoefficients
                catch ex
                    warning(warningIDBefore.state,'MATLAB:nearlySingularMatrix')
                    error(ex.message)
                end
                warning(warningIDBefore.state,'MATLAB:nearlySingularMatrix')
            
                SQ_CoVar = sum( (obj.getOutputData - obj.BasisFctCoefficients(1)*obj.BasisFct{1}(varargin{1},obj.getInputData)).^2 );
            else
                warning('For %s the function LeastSquareFunctionBasisFct will do nothing (no nonlinearparameter)',obj.allowed{1});
            end
        otherwise
            error('%s was not implemented',obj.BasisFctType)
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
