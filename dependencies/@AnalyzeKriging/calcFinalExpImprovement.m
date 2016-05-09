function [ExpectedImprovement] = calcFinalExpImprovement(obj,prediction,u)
% [ExpectedImprovement] = calcFinalExpImprovement(obj,prediction,u)
%
% Input:
% - prediction ... prediction made by Kriging at the point of interest [nPointsX2]
% - u ... variable which follows a standard normal distribution [nPointsX1]
%
% Output:
% - ExpectedImprovement ... expected improvement at the point of interest
%                           [nPointsX1]
%
% You can set: 
%
% - DegreeOfExpectedImprovement ... Exploration factor (integer >=0) the
%   bigger the more the will the solver follow an exploration strategy
%   (Prediciton error get a heigher weight).
%   % Extreme cases: 
%     1. DegreeOfExpectedImprovement=0 only the cumulative probability of
%     improvement is considered (Phi(I/sigma(x))).(Local Search)
%     2. DegreeOfExpectedImprovement>>1 Only the prediction variance is
%     considered to be important. Prediction variance depends mainly on the
%     distance to given data points -> Space filling procedure
%
% You can get: -
%
% Note: Further Details can be found in the PhD thesis of Matthias Schonlau
% "Computer Experiments and Global Optimization" (1998)
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.
   
% Initialization
ExpectedImprovement = 0;
if mod(obj.DegreeOfExpectedImprovement,1)~=0&&obj.DegreeOfExpectedImprovement<0
    error('DegreeOfExpectedImprovement has to be a positive integer for calculation of the single-objective expected improvement')
end

% Calculate Recursively factors t_k
t_k_vec = zeros(obj.DegreeOfExpectedImprovement+1,size(u,1));
t_k_vec(1,:) = Phi(u)';
t_k_vec(2,:) = -phi(u)';

% iDegreeOfExpectedImprovement actual start with
% iDegreeOfExpectedImprovement=0 but since indices in matlab start with
% 1, (+1) is added to iDegreeOfExpectedImprovement
for iDegreeOfExpectedImprovement=2:obj.DegreeOfExpectedImprovement
    factor = (iDegreeOfExpectedImprovement - 1);
    t_k_vec(iDegreeOfExpectedImprovement+1,:) = -phi(u)'.*u'.^factor + factor*t_k_vec(iDegreeOfExpectedImprovement+1 - 2,:);
end

for iDegreeOfExpectedImprovement=0:obj.DegreeOfExpectedImprovement

    ExpectedImprovement = ExpectedImprovement +...
                          (-1)^iDegreeOfExpectedImprovement.*...
                          nchoosek(obj.DegreeOfExpectedImprovement,iDegreeOfExpectedImprovement).*...
                          u.^(obj.DegreeOfExpectedImprovement-iDegreeOfExpectedImprovement).*...
                          t_k_vec(iDegreeOfExpectedImprovement+1,:)';
end

% Final Output
ExpectedImprovement = prediction(:,2).^obj.DegreeOfExpectedImprovement.*ExpectedImprovement;


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
