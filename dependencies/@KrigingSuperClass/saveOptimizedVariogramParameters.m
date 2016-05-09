function [] = saveOptimizedVariogramParameters(obj,optPara)
% [] = saveOptimizedVariogramParameters(optPara)
% This function is in general called after estimating the covariogram
% parameters. It sets the values given by optPara to the covariogram model
% parameters: theta,sigma, sigmaError, p depending on the definition of the
% covariogram model
%
%  You can set: -
%  
%  You can get: 
% - CovariogramModelParameters: ... assigne parameter set
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

    switch obj.CovariogramModelChoice
        case 1
            obj.theta       = optPara(1);
            obj.sigma       = optPara(2);
            obj.sigmaError  = 1e-10;
        case {2,5,7}
            obj.theta       = optPara(1);
            obj.sigma       = optPara(2);
            obj.sigmaError  = optPara(3);
        case 3
            nTotal = obj.nInputVar;
            obj.theta       = optPara(1:nTotal);
            obj.p           = optPara(nTotal+1:2*nTotal);
            obj.sigma       = optPara(2*nTotal+1);
            obj.sigmaError  = optPara(2*nTotal+2);
        case {4,6,8}
            nTotal = obj.nInputVar;
            obj.theta       = optPara(1:nTotal);
            obj.sigma       = optPara(nTotal+1);
            obj.sigmaError  = optPara(nTotal+2);
        otherwise
            error('Non-acceptable model choice. The parameter CovariogramModelChoice = %i is not allowed',obj.CovariogramModelChoice);
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
