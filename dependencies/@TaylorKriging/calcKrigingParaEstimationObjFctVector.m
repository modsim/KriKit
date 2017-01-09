function [quality] = calcKrigingParaEstimationObjFctVector(obj,para)
% Calculates the quality of the parameter estimation with respect to the
% Kriging quality measurement: log10(sum(a_l^2))
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

    % Create extended Covariance Matrix
        % Initialization
    if isempty(obj.CovariogramMatrix)
        obj.calcCovariogramMatrix
    end

    if length(para)~=obj.getnBasisFctParameters
        warning('Number of given parameters (%i) and the needed number of parameters (%i) differes from each other',length(para),obj.getnBasisFctParameters)
    end
    if obj.nBasisFctDerivative<(obj.getnBasisFctParameters+(obj.TaylorExpansionOrder>=2)*(obj.getnBasisFctParameters+obj.getnBasisFctParameters*(obj.getnBasisFctParameters-1)/2))
        error('Not enough derivatives are provided')
    end
    C = obj.CovariogramMatrix(1:obj.getnExperiments,1:obj.getnExperiments);
%     minC=min(min(abs(C)));
%     C = C/minC; % Normalization for reducing numerical problems
    
    F = zeros(obj.getnExperiments,1+obj.nBasisFctDerivative);
    
        % Basis Function
    F(:,1) = obj.getBasisFct{1}(para,obj.InputData);
        % Derivatives
    for iDerivative=1:obj.nBasisFctDerivative
        F(:,iDerivative+1) = obj.getBasisFctDerivative{iDerivative}(para,obj.InputData);
    end
    
    if isempty(obj.getInvCoVar)
        obj.InvCoVar = inv(C);
    end

%     obj.InvCoVar(1:5,1:5)
    saveCoeff = (F'*(obj.InvCoVar*F))\(F'*(obj.InvCoVar*(obj.OutputData)));
    quality = saveCoeff(2:end);
    
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
