function EHVI_3D=calcEHIV_3D_2016(paretoPoints,referencePoint,mu,sd)
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
% This function is inspired by the code of Michael Emmerich and Andre
% Deutz, LIACS, Leiden University, 2010 
% Based on the paper:
% M. Emmerich, A.H. Deutz, J.W. Klinkenberg: The computation of the
% expected improvement in dominated hypervolume of Pareto front
% approximations , LIACS TR-4-2008, Leiden University, The Netherlands
% http://www.liacs.nl/~emmerich/TR-ExI.pdf
%
% Copyright 2014-2015: Lars Freier, Eric von Lieres
% See the license note at the end of the file.


%% New Code Vectorized
% timeVec = tic;
% Make sure that only pareto optimal points are provided
paretoPoints = determineParetoSet_Mex(paretoPoints);
% Use last variable for sorting (easier for interpretation)
paretoPointsSort = sortrows(paretoPoints,3);
% Add Point. Needed for proper calculation

% paretoPointsSort = [-inf,-inf,referencePoint(3);paretoPointsSort;referencePoint(1),referencePoint(2),-inf];
paretoPointsSort = [-inf,-inf,paretoPointsSort(1,3);paretoPointsSort;referencePoint(1),referencePoint(2),referencePoint(2)];
% Sort with descend order
% paretoPointsSort = [referencePoint(1),referencePoint(2),-inf;paretoPointsSort;-inf,-inf,referencePoint(3)];
% paretoPointsSort = paretoPointsSort(end:-1:1,:);
nPoints = size(paretoPoints,1);

ParetoFront2D = referencePoint(1:2);
EHVI_3D = 0;
for iSection = 2:(nPoints+2)
    % Update 2D Pareto Curve
    if iSection>2
        ParetoFront2D = determineParetoSet_Mex(paretoPointsSort(2:iSection-1,1:2));
    end
    
    EHVI_2D = calcEHIV_2D(ParetoFront2D,referencePoint(1:2),mu(1:2),sd(1:2));
    
%     if iSection<(nPoints+2)
    if iSection<(nPoints+2)
        part1 = (paretoPointsSort(iSection,3)-paretoPointsSort(iSection-1,3)).*...
                 gausscdf((paretoPointsSort(iSection,3)-mu(3))/(sd(3))).*...
                 EHVI_2D;
    else
        part1 = 0;
    end
    
    if iSection>2
        part2 = (gaussEI(mu(3),sd(3),paretoPointsSort(iSection,3),paretoPointsSort(iSection,3))-...
                 gaussEI(mu(3),sd(3),paretoPointsSort(iSection,3),paretoPointsSort(iSection-1,3))).*...
                 EHVI_2D;
    else
        part2 =  gaussEI(mu(3),sd(3),paretoPointsSort(iSection,3),paretoPointsSort(iSection,3)).*...
                 EHVI_2D;
    end
    
    EHVI_3D = EHVI_3D + part1 + part2;
end


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