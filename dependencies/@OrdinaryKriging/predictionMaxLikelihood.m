function [output] = predictionMaxLikelihood(obj,input)
% [Output] = predictionMaxLikelihood(obj,input,significanceLevel)
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.
    % This Function Predicts the Outpus Values for the Given Input Points
    % input ... data points where the function value shall be predicted (nSamplePoints x nInputVar Matrix)
    % Output: Is a matrix nSamplePointsX2. nSamplePoints is the number of given sample points
    %         First column contains the actual prediction 
    %         Second column contains the standard deviation 
    
    
    
    % Calculate the between input and the measure sample points
    nInput = size(input,1);
    r = zeros(obj.nExperiments,nInput);
    for(iInput = 1 : nInput)
        for iExperiments = 1:obj.nExperiments
            % Calculate correlation between the one sample point which is to predict 
            % and the measured sample point.
            r(iExperiments,iInput) = obj.CovarModel(abs(input(iInput,:)-obj.InputData(iExperiments,:)));
%             distanceR(iExperiments,1:2) = (abs(input(iInput,:)-obj.InputData(iExperiments,:)));
        end
    end
    
    % Initialization
    output = zeros(nInput,2);
    
    % Calculate Prediction Values
    r2 = r'*obj.InvCoVar;
    output(:,1) = obj.MuMaxLikelihood + r2*(obj.OutputData-ones(obj.nExperiments,1)*obj.MuMaxLikelihood);
    
    % Calculate the standard deviation
%     r3 = r2*r;
    for(iInput = 1 : nInput)
        output(iInput,2) = (1-r2(iInput,:)*r(:,iInput)+((1-r2(iInput,:)*r(:,iInput)).^2 )./...
            (ones(1,obj.nExperiments)*obj.InvCoVar*ones(obj.nExperiments,1)));
    end
    output(:,2) = obj.SigmaMaxLikelihood*output(:,2);
    
    neg = find(output(:,2)<0);
    if ~isempty(neg)
        warning('%i negative standard deviations. Minimal SD is %d',length(neg),min(output(:,2)));
    end
    output(:,2) = abs(output(:,2));
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
