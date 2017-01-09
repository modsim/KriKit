function EHIV=calcEHIV_2D_2016(paretoPoints,referencePoint,mu,sd)
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


%% New Code
timeNonVec = tic;
% Make sure that only pareto optimal points are provided
paretoPoints = determineParetoSet_Mex(paretoPoints);
paretoPointsSort = sortrows(paretoPoints);
paretoPointsSort = [-inf,referencePoint(2);paretoPointsSort;referencePoint(1),-inf];
paretoPointsSort = paretoPointsSort(end:-1:1,:);

nSamples=size(paretoPoints,1);
nObj = size(paretoPoints,2);


EHIV = 0;
for iSection = 2:(nSamples+2)
    if iSection<(nSamples+2)
        part1 = (paretoPointsSort(iSection-1,1)-paretoPointsSort(iSection,1)).*...
                 gausscdf((paretoPointsSort(iSection)-mu(1))/(sd(1))).*...
                 gaussEI(mu(2),sd(2),paretoPointsSort(iSection,2),paretoPointsSort(iSection,2));
    else
        part1 = 0;
    end
         
    part2 = (gaussEI(mu(1),sd(1),paretoPointsSort(iSection-1,1),paretoPointsSort(iSection-1,1))-...
             gaussEI(mu(1),sd(1),paretoPointsSort(iSection-1,1),paretoPointsSort(iSection,1))).*...
            gaussEI(mu(2),sd(2),paretoPointsSort(iSection,2),paretoPointsSort(iSection,2));
        
    EHIV = EHIV + part1 + part2;
    
    fprintf('iSection: %i - x1(-1): %g - x1: %g - x2: %g - EHIV: %g - part1: %g - part2: %g\n',...
             iSection-1,paretoPointsSort(iSection-1,1),paretoPointsSort(iSection,1),paretoPointsSort(iSection,2),...
             EHIV,part1,part2)
    
end
fprintf('EHIV: %g - time: %g\n',EHIV,toc(timeNonVec))

%% New Code Vectorized
timeVec = tic;
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

fprintf('EHIV: %g - time: %g\n',EHIV,toc(timeVec))

%% Old Code
timeOld = tic;

% Make sure that only pareto optimal points are provided
paretoPoints = determineParetoSet_Mex(paretoPoints);
sortByEachDimension = sort(paretoPoints);

nSamples=size(paretoPoints,1);
nObj = size(paretoPoints,2);
nCellsPriori = (nSamples+1)^nObj;
nCells = (nSamples+1)*(nSamples+2)/2;

% For better clarity: sort every variable with increasing value
valuesCell = cell(nObj,1);
for iObj=1:nObj
    valuesCell{iObj} = sortByEachDimension(:,iObj);
end

% Define grid which represent the vertices of eahc considered rectangle
[i1,i2] = ndgrid(0:nSamples,0:nSamples);
indexMatrix =[i1(:),i2(:)];

% Samples in the recangle with these lower bounds to not lead to any
% contribution to the EIHV as anypoint in these areas are dominated by
% current pareto front
doNotUse = indexMatrix(:,1)>(nSamples-indexMatrix(:,2)); 
indexMatrix(doNotUse,:)=[];

% For more clarity, save lower and upper bounds of each reactangle
Cl_Matrix = zeros(nCells,nObj);
Cu_Matrix = zeros(nCells,nObj);
dominateSamplesCell = cell(nSamples-1,nObj);

for iIndex=0:nSamples
    for iObj=1:nObj
        if iIndex==0
            Cl_Matrix(indexMatrix(:,iObj)==0,iObj) = -inf;
        else
            Cl_Matrix(indexMatrix(:,iObj)==iIndex,iObj) = valuesCell{iObj}(iIndex);
        end
    end
end

for iIndex=0:nSamples
    for iObj=1:nObj
        if iIndex==nSamples
            Cu_Matrix(indexMatrix(:,iObj)==iIndex,iObj) = referencePoint(iObj);
        else
            Cu_Matrix(indexMatrix(:,iObj)==iIndex,iObj) = valuesCell{iObj}(iIndex+1);
        end
    end
end

% Also save values which are dominated by each sample points in each
% direction
for iIndex=0:nSamples-1
    for iObj=1:nObj
        dominateSamplesCell{iIndex+1,iObj} = paretoPoints(:,iObj)>=valuesCell{iObj}(iIndex+1);
    end
end

% Find reference points fmaxMatrix for each considered rectangle
fmaxMatrix = bsxfun(@times,ones(nCells,nObj),referencePoint);
for iCell=1:nCells
    for iObj=1:nObj
        if iObj==1
            indexChosenParetoPoint = paretoPoints(:,2)==Cl_Matrix(iCell,2);
        else
            indexChosenParetoPoint = paretoPoints(:,1)==Cl_Matrix(iCell,1);
        end
        if ~any(indexChosenParetoPoint)
            fmaxMatrix(iCell,iObj) = referencePoint(iObj);
        else
            fmaxMatrix(iCell,iObj) = min(paretoPoints(indexChosenParetoPoint,iObj));
        end
    end
end

% sampleIsDominated = false(nSamples,1);
sPlus = zeros(nCells,1);


for iCell=1:nCells
        % Identify sample points which are dominated by upper corner of
        % current cell
        if indexMatrix(iCell,1)==nSamples||indexMatrix(iCell,2)==nSamples
            sampleIsDominated = false(nSamples,1);
        else
            sampleIsDominated = dominateSamplesCell{indexMatrix(iCell,1)+1,1}&dominateSamplesCell{indexMatrix(iCell,2)+1,2};
        end
        

        if ~any(sampleIsDominated)
            sPlus(iCell)=0;
        else
            sPlus(iCell) = Hypervolume_MEX(paretoPoints(sampleIsDominated,:),fmaxMatrix(iCell,:));
        end
end

%Marginal integration 
Psi = marginalEI(fmaxMatrix,Cu_Matrix,mu,sd) - marginalEI(fmaxMatrix,Cl_Matrix,mu,sd);

%Cumulative Gaussian over length for correction constant
Cu_Norm = bsxfun(@rdivide,bsxfun(@minus,Cu_Matrix,mu),sd);
Cl_Norm = bsxfun(@rdivide,bsxfun(@minus,Cl_Matrix,mu),sd);

GaussCDF = gausscdf(Cu_Norm) - gausscdf(Cl_Norm);
%ExI Kontribution fuer die aktuelle Zelle
cellEI = prod(Psi,2)-sPlus.*prod(GaussCDF,2);


EHIV=sum(sum(cellEI));

fprintf('EHIV: %g - time: %g\n',EHIV,toc(timeOld))

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