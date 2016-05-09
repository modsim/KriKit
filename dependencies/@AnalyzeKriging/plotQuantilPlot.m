function [PearsonCoeff,SumOfSquares] = plotQuantilPlot(obj,varargin)
% [] = plotQuantilePlot(KrigingObjectIndex) 
%
% This function creates a quantile plot for a better evaluation of the
% Kriging estimation results. Data points should follow the black linear
% curve. 
%
% Background: 
% The difference between the Kriging prediction and the actual mesurements
% (residual) should be unbiased with minimal variance. In the framework of
% Kriging, it is assumed that the measurment noisy is identical for all
% data points.(something between sigma_Error and sigma_Error+sigma). 
% Consequently the centralized residual weighted by the standard deviation
% should follow a standard normal distribution
%
% Input:
% KrigingObjectIndex ... Index of Kriging object of interest
%
% Output:
% PearsonCoeff ... Pearson coefficient as quantifier for linear correlation
%                  be x- and y-axis
% SumOfSquares ... Sum of squared residuals (=sum(((x-value)-(y-value))^2) )
%
% You can set: -
%
% You can get: -
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

%% Presteps
KrigingObjectIndex = varargin{1};

% Do Prediction
predictionOutput = obj.KrigingObjects{KrigingObjectIndex}.prediction(...
                        obj.KrigingObjects{KrigingObjectIndex}.getInputData);
                    
% Calculate coummulative distribution
nExperiments = obj.KrigingObjects{KrigingObjectIndex}.getnExperiments;
residuals = (predictionOutput(:,1)-obj.KrigingObjects{KrigingObjectIndex}.getOutputData)./predictionOutput(:,2);

residualsSort = sort(residuals);
cumProb = zeros(nExperiments,1);

% Calculate for each risidual the cum. prob.
for iRisdual = 1:nExperiments
    cumProb(iRisdual) = sum(residuals<=residualsSort(iRisdual))/(nExperiments);
end

%% Calculate the ideal quantiles
fQuantile = @(q)sqrt(2)*erfinv(2*q-1);

%% Create Plot
% If process is ideal, than the centralized residual should be equal to
% the quantile of the standard normal distribution given the actual
% cumulative probability
figure()
hold on
plot(fQuantile(cumProb(1:end-1)),residualsSort(1:end-1),'ko','MarkerFaceColor',[255/255 102/255 0/255])

% Straight line
plot(fQuantile(cumProb(1:end-1)),fQuantile(cumProb(1:end-1)),'k','LineWidth',3)
xlabel('Ideal Cumulative Probability','FontSize',20)
ylabel('Actual Cumulative Probability','FontSize',20)

%% Get Output
pearsonCoeffMatrix = corrcoef(fQuantile(cumProb(1:end-1)),residualsSort(1:end-1));
PearsonCoeff = pearsonCoeffMatrix(1,2);
SumOfSquares = sum((fQuantile(cumProb(1:end-1)) - residualsSort(1:end-1)).^2);
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
