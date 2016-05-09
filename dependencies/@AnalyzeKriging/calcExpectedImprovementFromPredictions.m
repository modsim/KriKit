function [ExpectedImprovement]=calcExpectedImprovementFromPredictions(obj,varargin)
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.
% [ExpectedImprovement]=calcExpectedImprovement(obj,KrigingObjectIndex,predictions)
% KrigingObjectIndex ... indicates the which Kriging Objective
% is considered
% prediction ... nPoints X 1 matrix containing Kriging predictions

    KrigingObjectIndex = varargin{1};
    if length(KrigingObjectIndex)~=1
        error('KrigingObjectIndex must be a scalar')
    end
    
    % Make Prediction
    prediction = varargin{2};
    prediction(:,1) = -obj.MinMax(1)*prediction(:,1);

    % Calculate the expected improvements at the particular sample
    % locations
        % Best Measurement so far
    %             minOutput = min(obj.MinMax(1)*obj.KrigingObjects{KrigingObjectIndex}.getOutputData);
     outputData = -obj.MinMax(KrigingObjectIndex)*obj.KrigingObjects{KrigingObjectIndex}.getOutputData;
    [~,indexMin] = min(outputData);
    bestOutput = outputData(indexMin);
        % Argument for cumulative distribution function
    u = (bestOutput - prediction(:,1))./prediction(:,2);
        % Calculate expected improvement
            % if statistic toolbox is installed, the following command can be used:
            % ExpectedImprovement = (minOutput - prediction(:,1)).*normcdf(u,0,1) + prediction(:,2).*normpdf(u,0,1);
            % ExpectedImprovement = (minOutput - prediction(:,1)).*(1/2+1/2*erf(u/sqrt(2)))+...
            %                         prediction(:,2)./sqrt(2*pi).*exp(-u.^2/2);
   switch obj.DegreeOfExpectedImprovement
        case 0
            ExpectedImprovement = Phi(u);
        case 1
            ExpectedImprovement = prediction(:,2).*(u.*Phi(u)+phi(u));
        case 2
            ExpectedImprovement = prediction(:,2).^2.*((u.^2+1).*Phi(u)  +u.*phi(u));
        case 3
            ExpectedImprovement = prediction(:,2).^3.*((u.^3+3*u).*Phi(u) + (2+u.^2).*phi(u));
       otherwise
            [ExpectedImprovement] = calcFinalExpImprovement(obj,prediction,u);
    end

    %% Nested Function
    function [cumProb]=Phi(x)
        cumProb = (1/2+1/2*erf(x/sqrt(2)));
    end
    function [prob]=phi(x)
        prob = 1./sqrt(2*pi).*exp(-x.^2/2);
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
