function [ ] = makePrework(obj)
% [] = makePrework()
% This function calculates the Covariogram parameters
%
% You can set: 
% - UseMatlabRegressionGP ... if true, Mathlab statistical Toolbox is used
%                             for parameter erstimation
% - UseSolver ... decide between different solvers (only if "UseMatlabRegressionGP"==false): 
%      1 ... fmincon for further information see documentation of
%      "solveLeastSquareCovariogram()"
%      2 ... Genetic Algorithm by Andrew Chipperfield, for further
%      information see documentation of "solveLeastSquareCovariogramGA2()"
%      3 ... Genetic Algorithm included in the Global Optimization Toolbox
%      of Matlab. For information see documentation of
%      "solveLeastSquareCovariogramGA()" 
%
% 
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

%     switch obj.MakePreworkType
%         case {1,3,4}
    if obj.UseMatlabRegressionGP
        obj.generateRegressionGPModel;
    else
        % Calculate the experimental Semi-Variogramm
        if isempty(obj.getDistInput)
            obj.estimateVariance();
        end

        % Find Approximation of the Semi-Variogramm
        switch obj.UseSolver
            case 1
                obj.solveLeastSquareCovariogram();
            case 2
                obj.solveLeastSquareCovariogramGA2();
            case 3
                obj.solveLeastSquareCovariogramGA();
            otherwise
                error('Flag UseSolver "%d" was not defined',obj.UseSolver)
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
