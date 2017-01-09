function [negLog] = MarginalLikelihoodCovarPara(obj,varargin)
% [negLog] = MarginalLikelihoodCovarPara(parameters)
% Calculates the marginal likelihood
% p(output|input,covParam,covModel)=int(p(output|input,basisFunctionParmeter,covModel)*p(basisFunctionParmeter|covParam,covModel))
%
% negLog=-log(p(output|input,covParam,covModel)) = -1/2*(-(output - basisFct)'*(C\(output-basisFct)) -log(det(covMatrix))-nOutput*log(2*pi));
%
% "It is primarily the marginal likelihood involving the integral
% over the parameter space which distinguishes the Bayesian scheme of inference
% from other schemes based on optimization. It is a property of the marginal
% likelihood that it automatically incorporates a trade-off between model fit and
% model complexity. This is the reason why the marginal likelihood is valuable
% in solving the model selection problem." (Rasmussen 2006: "Gaussian
% Processes for Machine Learning")
%
% You can set: 
% - CovariogramModelChoice ... decied for covariogram model
%
% You can get: -
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

% Define the model parameter values 
saveOptimizedVariogramParameters(obj,varargin{1});
calcCovariogramMatrix(obj);
estimateBasisFctCoefficients(obj)
C = obj.CovariogramMatrix(1:obj.nExperiments,1:obj.nExperiments);
detC = det(C);
numericalLimit = realmin;
% if detC<=numericalLimit
if isinf(detC)||detC==0
    detC = numericalLimit;
end

lC = lu(C);

[basisFctMatrix] = estimatedBasisFunction();
basisFctEval = basisFctMatrix'*obj.getBasisFctCoefficients;
vecForMulti = (obj.OutputData-basisFctEval);
negLog = -(-1/2*vecForMulti'*(C\vecForMulti) - 1/2*log(detC)-obj.nExperiments/2*log(2*pi));
% negLog = -(-1/2*vecForMulti'*(C\vecForMulti) - sqrt(obj.nExperiments)*log(min(abs(diag(lC))))-obj.nExperiments/2*log(2*pi));

%% Nested Function
% -------------------------------------------------------------------------
function [basisFctMatrix] = estimatedBasisFunction()
    % Extend Matrix by Basis Functions
    basisFctMatrix = zeros(obj.nBasisFct,obj.nExperiments);
    for iBasis = 1 : obj.nBasisFct
        basisFctMatrix(iBasis,:) = (obj.BasisFct{iBasis}(obj.BasisFctParameters,obj.getInputData))';
    end
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
