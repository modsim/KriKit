function [Output] = prediction(obj,input)
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.
    % [Output] = prediction(obj,input,significanceLevel)
    % This Function Predicts the Outpus Values for the Given Input Points
    % input ... data points where the function value shall be predicted (nxInputVar Matrix)
    % significanceLevel ... probability of false negative estimation
    % Output: Is a matrix nX3. n is the number of given data points
    %         First column contains the lowerbounds of the estimation
    %         Second column contains the actual prediction 
    %         Third column contains the upper bound of the estimation

    if (obj.checkInvVariogram==0)
        obj.calcInverseVariogram();
        obj.checkInvVariogram = 1;
        disp('Inverse was calculated');
       
    end
    
    VarGramm = zeros(obj.nExperiments+1,length(input));
% ------------------------------------------------------------------------    
    % Calculate Variogram Values the Given Inputs
    for(iData = 1 : obj.nExperiments)
        for(iInputs=1:size(input,1))
            VarGramm(iData,iInputs) = obj.CovarModel(0,1) - obj.CovarModel(norm(input(iInputs,:)-obj.InputData(iData,:),2),1);
%             VarGramm(nData,nInputs) = theta - obj.CovarModel(norm(input(nInputs,:)-obj.InputData(nData,:),2));
        end
    end
    
    
    % To allow multiple measurements at one point
%     VarGramm = VarGramm - obj.sigmaError;
    VarGramm(obj.nExperiments+1,:) = ones(1,size(input,1));
% ------------------------------------------------------------------------
    % Calculate Actual Prediction
    weights = obj.InvVariogram*VarGramm;
    estimation = [obj.OutputData;0]'*weights;
    estimation = estimation';
% ------------------------------------------------------------------------
    % Calculate the Variance at each Point
    sigmaEstimation = VarGramm'*weights;
        % Only the diagonal are the number we are looking for
        % var'*InvVariogram*var
        % var .. vector of variogram between new data points and existing
        % ones
        % InvVariogram ... inverse of the extended variogramm matrix
    sigmaEstimation = diag(sigmaEstimation);
    sigmaEstimation = (sigmaEstimation).^(1/2);
    % Check for imaginary numbers
    testImag = find(imag(sigmaEstimation)~=0);
    if(~isempty(testImag))
        warning('%i number of the standard deviations are imaginary',length(testImag))
        sigmaEstimation = real(sigmaEstimation);
    end
% ------------------------------------------------------------------------
    % Calculate Quantile: (1+certainityLevel)/2
%     try
% %         keyboard
%         quantileValue=norminv((1+1-significanceLevel)/2);
%     catch
%         warning('It seems that Statistic ToolBox is not installed. Significance level of 0.05 is assumed.\n')
%         quantileValue = 1.96;
%     end
    
    Output = [estimation,sigmaEstimation];
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
