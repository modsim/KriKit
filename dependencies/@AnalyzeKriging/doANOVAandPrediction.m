function [] = doANOVAandPrediction(obj,KrigingObjectIndex,X,inputMatrix,df)
% This Function does the main work of the ANOVA
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

%% Initialization
nInputVar = obj.KrigingObjects{KrigingObjectIndex}.getnInputVar;
nData = size(inputMatrix,1);

%% Calculate Kriging Interpolation
useGPRBackup = false;
if obj.KrigingObjects{KrigingObjectIndex}.getUseMatlabRegressionGP
    useGPRBackup = true;
    obj.KrigingObjects{KrigingObjectIndex}.setUseMatlabRegressionGP(false);
end
krigingPrediction = obj.KrigingObjects{KrigingObjectIndex}.prediction(inputMatrix);
obj.KrigingObjects{KrigingObjectIndex}.setUseMatlabRegressionGP(useGPRBackup);

krigingMeanValue = krigingPrediction(:,1);
krigingSigma = krigingPrediction(:,2);

%% Do Analysis
obj.KrigingObjects{KrigingObjectIndex}.calcCovarMatrixOfPredictions(inputMatrix)
% sigmaMatrix = (bsxfun(@times,eye(nData),krigingSigma)).^2;
sigmaMatrix = obj.KrigingObjects{KrigingObjectIndex}.getCovarMatrixOfPredictions;


% Regression
obj.ANOVACoefficients = (X'*(sigmaMatrix\X))\(X'*(sigmaMatrix\krigingMeanValue));
% beta = (X'*(X))\(X'*(krigingMeanValue));
% varExp = sum((krigingMeanValue - X*beta).^2)/(nData-nPara);
% covBeta = varExp*inv(X'*X);
obj.ANOVACovMatrix = inv(X'*(sigmaMatrix\X));
obj.ANOVAStdOfCoefficients = real(sqrt(diag(obj.ANOVACovMatrix)));

% T-Test
obj.ANOVATvalue=obj.ANOVACoefficients./obj.ANOVAStdOfCoefficients;

obj.ANOVAPvalue = 1- (1 - betainc(df./(df+obj.ANOVATvalue.^2),0.5*df,0.5));

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
