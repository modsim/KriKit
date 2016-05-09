function [] = InitializeCovariogramParameters(obj)
% [] = InitializeCovariogramParameters()
% This function initialize the covariogram parameters if no initial
% parameter set is yet defined (InitialCovariogramParameters is empty)
% This is need for functions such as "solveLeastSquareCovariogram"
%
% You can set: 
% - CovariogramModelChoice ... decide for the covariogram model
%
% You can get: 
% - InitialCovariogramParameters ... initialization set
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

    if isempty(obj.InitialCovariogramParameters)
        switch obj.CovariogramModelChoice
            case {1}
                obj.InitialCovariogramParameters = [1 1];
            case {2,5,7}
                obj.InitialCovariogramParameters = [1 1 1e-10];
            case 3
                nTotal = obj.nInputVar;
                obj.InitialCovariogramParameters = [ones(1,nTotal) ones(1,nTotal)*1 1e-10 1e-10];
            case {4,6,8}
                nTotal = obj.nInputVar;
                obj.InitialCovariogramParameters = [ones(1,nTotal) 1e-10 1e-10];
            otherwise
                error('Non-acceptable model choice. The parameter CovariogramModelChoice = %i is not allowed',obj.CovariogramModelChoice);
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
