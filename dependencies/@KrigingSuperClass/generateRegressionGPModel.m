function [] = generateRegressionGPModel(obj)
% Generate a Gaussian Process Regression model for further Kriging estimation
% This command replaces "makePrework" when "UseMatlabRegressionGP=true"
%
% You can set:
% - CovariogramModelChoice ... decide for the used covariogram model. Sofar
%   only the following model can be used for GPR: 'squaredExpWithoutNugget',
%   'squaredExponential', 'squaredExponentialComplex',
%   'ardsquaredexponential', 'matern32', 'ardmatern32', 'matern52',
%   'ardmatern52'
% - BasisFct ...represent the detemrinic trend for more information see
%               documentation of "setBasisFct()". Currently only constant
%               or linear basis function can be used
%
% % You can get:
% - GPR_Model ... Resulting Gaussian process regression model
% - CovariogramModelParameters ... p,theta, sigma, sigmaError
%
% Note: For using this function Matlab2015b or newer is ncessary and the
% stastitic toolbox has to be avaiable
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

if sum([2,4:8]==obj.getCovariogramModelChoice)~=1
    error('For Matlab Gaussian Process Regression Modeling only CovariogramModelChoice 2&4:8 are allowed')
end

% Define Basis Function 
warningMessage = 'generateRegressionGPModel only supports polynomials with degree of 0 and 1. Degree of 0 is chosen by default';
typeGPR_BasisFct = 'constant';
if strcmp(obj.BasisFctType,'polynomial')
    switch obj.MaxDegree
        case 0
            typeGPR_BasisFct = 'constant';
        case 1
            typeGPR_BasisFct = 'linear';
        otherwise
            warning(warningMessage)
    end
else
    warning(warningMessage)
end

models = {'squaredExpWithoutNugget','squaredExponential','squaredExponentialComplex','ardsquaredexponential','matern32','ardmatern32','matern52','ardmatern52'};

% sUse Data As they are
obj.GPR_Model = fitrgp(obj.InputData,obj.OutputData,'KernelFunction',models{obj.getCovariogramModelChoice},'BasisFunction',typeGPR_BasisFct);

p = obj.GPR_Model.KernelInformation.KernelParameters;
sigma = obj.GPR_Model.Sigma;
obj.setCovariogramModelParameters([p',sigma])

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
