function [optLambda]=calcLambdas(obj)
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

    % Initial lambdas
    Cinv = calcCinv;
    para = obj.getBasisFctParameters;
    k=obj.getnBasisFctParameters;
    n=obj.getnExperiments;
    allComb = 1+k*2+k*(k-1)/2;
    iniLambda = zeros(allComb*n+allComb*k+k*k,1);
    iniLambda(1:n*(obj.nBasisFctDerivative+1))=reshape(Cinv(1:n,n+1:end),n*(obj.nBasisFctDerivative+1),1);
    options = obj.setOptionsFminCon;
    optLambda=fmincon(@obj.calcLamdaObjFct,iniLambda,[],[],[],[],[],[],@obj.calcLamdaContraints,options);

%% Nested Functions
    function [Cinv]=calcCinv()
        para = obj.getBasisFctParameters;
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
        C = zeros(obj.getnExperiments+1+obj.nBasisFctDerivative);
        C(1:obj.getnExperiments,1:obj.getnExperiments) = ...
                            obj.CovariogramMatrix(1:obj.getnExperiments,1:obj.getnExperiments);

            % Basis Function
        C(1:obj.getnExperiments,obj.getnExperiments+1) = obj.getBasisFct{1}(para,obj.getInputData);
            % Derivatives
        for iDerivative=1:obj.nBasisFctDerivative
            C(1:obj.getnExperiments,obj.getnExperiments+iDerivative+1) = obj.getBasisFctDerivative{iDerivative}(para,obj.getInputData);
        end

        C=triu(C)+triu(C)'-bsxfun(@times,eye(length(C)),diag(C));
        Cinv = inv(C);
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
