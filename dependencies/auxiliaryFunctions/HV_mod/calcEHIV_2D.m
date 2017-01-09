function EHIV=calcEHIV_2D(paretoPoints,referencePoint,mu,sd)
% EHIV=calcEHIV_2D(paretoPoints,referencePoint,mu,sd)
% 
% This function calculate the expected hypervolume improvement of a Kriging
% model based on data point with the Pareto front given in paretoPoints.
%
% Input:
% - paretoPoints ... contains the values of the Pareto front
%                    (nParetoPointsX2) 
% - referencePoint ... Reference point which is used for the hyper volume
%                      calculation. (1X2)
% - mu ... contains the Kriging prediction of the point of interest (1X2)
% - sd ... contains the Kriging prediction error (1X2)
%
% Output:
% - EHIV ... expected hypervolume improvement
%
% This function is inspired by the publication of Michael Emmerich et Al:
%
% M. Emmerich, K. Yang, A. Deutz, H. Wang, and C. M. Fonseca, “A
% multicriteria generalization of Bayesian global optimization,” in
% Advances in Stochastic and Deterministic Global Optimization (P. M.
% Pardalos, A. Zhigljavsky, and J. Žilinskas, eds.), vol. 107 of Springer
% Optimization and Its Applications, Springer International Publishing,
% 2016. In press.
%
% Copyright 2014-2015: Lars Freier, Eric von Lieres
% See the license note at the end of the file.


%% New Code Vectorized
% timeVec = tic;
% Make sure that only pareto optimal points are provided
paretoPoints = determineParetoSet_Mex(paretoPoints);
paretoPointsSort = sortrows(paretoPoints);
paretoPointsSort = [-inf,referencePoint(2);paretoPointsSort;referencePoint(1),-inf];
paretoPointsSort = paretoPointsSort(end:-1:1,:);

part1 = (paretoPointsSort(1:end-2,1)-paretoPointsSort(2:end-1,1)).*...
         gausscdf((paretoPointsSort(2:end-1,1)-mu(1))/(sd(1))).*...
         gaussEI(mu(2),sd(2),paretoPointsSort(2:end-1,2),paretoPointsSort(2:end-1,2));
part1 = [part1;0];

part2 = (gaussEI(mu(1),sd(1),paretoPointsSort(1:end-1,1),paretoPointsSort(1:end-1,1))-...
             gaussEI(mu(1),sd(1),paretoPointsSort(1:end-1,1),paretoPointsSort(2:end,1))).*...
            gaussEI(mu(2),sd(2),paretoPointsSort(2:end,2),paretoPointsSort(2:end,2));
EHIV = sum(part1+part2);

% fprintf('EHIV: %g - time: %g\n',EHIV,toc(timeVec))
% =============================================================================
%  KriKit - Kriging toolKit
%  
%  Copyright 2014-2015: Lars Freier(1), Eric von Lieres(1)
%                                      
%    (1)Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
%  
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================