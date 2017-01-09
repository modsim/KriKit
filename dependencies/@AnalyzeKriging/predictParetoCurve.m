function [expectedParetoCurve,deviationParetoCurve,pCover,gridOutput] = predictParetoCurve(obj,varargin)
% [expectedParetoCurve,deviationParetoCurve,pCover,gridOutput] = predictParetoCurve(krigingObjVector,nRealizations,nSampleLocations,nGridPointEachDirection,inEqualityCellArray,OutputConstraints)
%
% This function estimated the expected pareto curve pareto based on
% conditional simulation results. Theory is described in detail in is
% Binois2016: "Quantifying uncertainty on Pareto fronts with Gaussian
% process conditional simulations"
%
% Short Desciption: 
% nSampleLocations samples are distributed over the input space, using
% lattin hyper cuve design. Conditional simulations are performed at these
% location for all kriging object defined in krigingObjVector. For further
% details about conditional simulation see documentation of
% "doConditionalSimulation_ResAlgInEq()". 
%
% NOTE: If no licence for the statistical tool box is available, samples
% are distributed over a nSampleLocations-grid. Inequalities can be
% provided via "inEqualityCellArray". Where should be a handle to a
% function of the form
% [boolIsValidVec]=inEqualityCellArray(inputDataMatrix). Where
% inputDataMatrix (nInputDataXnInputVar) is a matrix containing the test
% points and boolIsValidVec is a vector where an entry at poistion i of
% true indicates that the input point at row i of inputDataMatrix holds the
% inequality condition. All points which do not hold the inequlity
% conditions are removed from the grid.
%
% A nGridPointEachDirection^nKrigingObjects grid is afterwards created over
% the output space. For each output grid point, the frequency, pCover, is
% calculated representing how often the point is dominated by a realization
% r created by the conditional simulation. Grid points with same value of
% pCover are called the pCover-Quantile. The expected Pareto curve contains
% output grid point of the mu-qunatile which satisfies HV(mu-quntile) =
% mean(HV(r)), where mean(HV(r)) represents the mean value of all
% Hypervolume spanned by the conditioned simulations and HV(mu-quntile) the
% Hypervolume spanned by the mu-quntile. (Vorob'ev expectation)
%
% Furthermore, Vorob'ev deviation is calculated in this function. It
% describes the frequency that a point is dominated by a realization and
% not by the expected pareto curve or vice versa. It might also be
% interpreted as the probability of point to be part of the pareto curve.
%
% Input: 
% - krigingObjVector ... indices of the kriging objects which are of
%                        interest
% - nRealizations ... number of conditional simulations for each kriging
%                     object
% - nSampleLocations ... number of samples which are distributed over the
%                        input space. If a nSampleLocations-grid is
%                        used, each input variable has
%                        nSampleLocations^/1/nInputVar) different levels
% - nGridPointEachDirection ... number of level for each kriging objective
%                               which is used for creating the
%                               nGridPointEachDirection^nKrigingObjects
%                               (output)-grid
% - inEqualityCellArray ... handle refering to a function which takes as
%                           input a (nSamplePointsXnInputVar)-Matrix at put
%                           out a (nSamplePointsX1)-boolean vector. Where
%                           each entry indicates if the associated sample
%                           point is valid (true) or not (false).
% - OutputConstraints ... a matrix (nKrigingObjectsX2) containing the
%                         feasible lower and upperbound for the output
%                         space which shall be investigated
% Output: 
% - expectedParetoCurve ... objective matrix containing the Pareto points
%                           (nParetoPointsXnKrigingObjects)
% - deviationParetoCurve ... probability Vorob'ev deviation over "grid"
%                            (nGridPointEachDirection^nKrigingObjectsXnKrigingObjects) 
% - gridOutput ... matrix containing the
%                  (output)-grid
%                  (nGridPointEachDirection^nKrigingObjectsXnKrigingObjects) 
%
% You can set:
% - setShowData ... Scan solution space also w.r.t. to the given output
%   data. This might needed if you want to plot the expected Pareto Curve
%   together with the data. Caution: This might lead to a ncombinatorical
%   explosion!
% - LB(UB)InputVarInterpolation ... restict the input space for the
%                                   conditional simulation
%  - RepeatDesign ... If "false" than a different lattin hypercube design
%                     is used for each realization (more roust but more
%                     expensive)
%
% You can get: - 
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

%% Initialization
krigingObjVector = varargin{1};
nObjects = length(krigingObjVector);
nRealizations = varargin{2};
% Assume equal number of input variables for each objective
nInputVar = obj.KrigingObjects{krigingObjVector(1)}.getnInputVar;
nSampleLocations = varargin{3};
if length(varargin)>3&&~isempty(varargin{4})
    nGridPointEachDirection = varargin{4};
    if length(varargin)>4
        inequalityConstraintInput = varargin{5};
    else
        inequalityConstraintInput = {};
    end
else
    nGridPointEachDirection = 1e6^(1/nObjects);
    inequalityConstraintInput = {};
end
if length(varargin)>4&&~isempty(varargin{5})
    inEqualityCellArray = varargin{5};
else
    inEqualityCellArray = {};
end
if length(varargin)>5&&~isempty(varargin{6})
    OutputConstraints = varargin{6};
else
    OutputConstraints = ones(nObjects,2);
    OutputConstraints(:,1) = -inf;
    OutputConstraints(:,2) = inf;
end


% If inequalityConstraintInput is given, nSampleLocations might change
defineBoundOfInputVar(obj,krigingObjVector(krigingObjVector(1)));
sampleGrid = createGrid;

%% #### Do Conditional Simulation #### 
realizationCurveMaxtrix = zeros(nSampleLocations,nObjects*nRealizations);
returnSamplePointsCell = zeros(nSampleLocations,nObjects*nRealizations*nInputVar);
for iKrigingObject = 1:nObjects
    if any(obj.getShowDetails)
        fprintf('Conditional Simulation for KrigingObject %i (%i of %i)\n',krigingObjVector(iKrigingObject),iKrigingObject,nObjects);
    end

    if ~isempty(inEqualityCellArray)
        [realizationCurve,returnSamplePoints] =...
                obj.KrigingObjects{krigingObjVector(iKrigingObject)}.doConditionalSimulation_ResAlgInEq(sampleGrid,inEqualityCellArray{iKrigingObject});
    else
        [realizationCurve,returnSamplePoints] =...
                obj.KrigingObjects{krigingObjVector(iKrigingObject)}.doConditionalSimulation_ResAlgInEq(sampleGrid);
    end

    columnIndicesRealization = (iKrigingObject-1)*nRealizations+1:iKrigingObject*nRealizations;
    columnIndicesSamples = (iKrigingObject-1)*nRealizations*nInputVar+1:iKrigingObject*nRealizations*nInputVar;
    realizationCurveMaxtrix(:,columnIndicesRealization) = realizationCurve;
    returnSamplePointsCell(:,columnIndicesSamples) = returnSamplePoints;
end

%% #### Find RNP (random non-dominated points) for each realization #### 
% Memory allociation
indexSet = zeros(1,nObjects);
paretoSets = cell(nRealizations,1);
nMemberPratoCurve = zeros(nRealizations,1);
totalParetoSet = [];

if any(obj.getShowDetails)
    fprintf('Collect Random Non-dominated Points (RNP)\n');
end

for iRealization = 1 : nRealizations
    for iObjChoice = 1:nObjects
        indexSet(iObjChoice) = (iObjChoice-1)*nRealizations + iRealization;
    end

    paretoSets{iRealization} = determineParetoSet_Mex(...
            bsxfun(@times,realizationCurveMaxtrix(:,indexSet),-obj.MinMax(krigingObjVector)) );
        
    % Save Pareto Set in Matrix
    totalParetoSet = [totalParetoSet;paretoSets{iRealization}];
    nMemberPratoCurve(iRealization) = size(paretoSets{iRealization},1);
end

%% #### Calculate pCover #### 
if any(obj.getShowDetails)
    fprintf('Calculate Cover Probability\n');
end

% Determine Range for Grid in the output space
[rangeToCheck] = determineRangeForGrid;

% Create Grid
if obj.ShowData
    gridX = createNDGRID(rangeToCheck(:,1),rangeToCheck(:,2),nGridPointEachDirection,totalParetoSet);
else
    gridX = createNDGRID(rangeToCheck(:,1),rangeToCheck(:,2),nGridPointEachDirection);
end

nGridPoints = size(gridX,1);

% Actual Calculation
pCover = mainMexCalcPcover(nGridPoints,...
              nRealizations,...
              nMemberPratoCurve,...
              nObjects,...
              totalParetoSet,...
              gridX,true);

%% #### Determine expected Pareto Curve #### 
if any(obj.getShowDetails)
    fprintf('Determine expected Pareto Curve\n');
end

% Calculate mean hyper volume
HV_Vector = zeros(nRealizations,1);
for iRealization = 1 : nRealizations
    if isempty(paretoSets{iRealization})
        error('paretoSets is empty. Check "MinMax" and "ReferencePointHyperVolume"')
    end
    HV_Vector(iRealization) = Hypervolume_MEX(paretoSets{iRealization},obj.ReferencePointHyperVolume(krigingObjVector).*(-obj.MinMax(krigingObjVector)) );
end
HV_Mean = mean(HV_Vector);


% Actual Calculation
[expectedParetoCurve] = determineExpectedParetoCurve;

%% #### Calculate Uncertainity (total) #### 
if any(obj.getShowDetails)
    fprintf('Calculate Uncertainity\n');
end
globalUncertainityVec =  zeros(nRealizations,1);
globalUncertainityVecNorm =  zeros(nRealizations,1);
for iRealization=1:nRealizations
    unionPoints = [expectedParetoCurve;paretoSets{iRealization}];
    HV_Union = Hypervolume_MEX(unionPoints,obj.ReferencePointHyperVolume(krigingObjVector).*(-obj.MinMax(krigingObjVector)));
    globalUncertainityVec(iRealization) = 2*HV_Union - HV_Quantile - HV_Vector(iRealization);
    globalUncertainityVecNorm(iRealization) = (2*HV_Union - HV_Quantile - HV_Vector(iRealization))/(2*HV_Union);
end
obj.GlobalParetoUncertainity = mean(globalUncertainityVec);
obj.GlobalParetoUncertainityNorm = mean(globalUncertainityVecNorm);

%% #### Calculate Uncertainity (single Point) #### 
deviationParetoCurve = mainMexCalcPcover(nGridPoints,...
              nRealizations,...
              nMemberPratoCurve,...
              nObjects,...
              totalParetoSet,...
              gridX,...
              false,...
              expectedParetoCurve,...
              pCover);
 deviationParetoCurve(1:end) = deviationParetoCurve(end:-1:1);
 pCover(1:end) = pCover(end:-1:1);
 
%% #### Convert everything back to the original space #### 
expectedParetoCurve = bsxfun(@times,expectedParetoCurve,-obj.MinMax(krigingObjVector));
gridOutput = bsxfun(@times,gridX,-obj.MinMax(krigingObjVector));

%% ------------ Nested Functions ------------ 
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
function [rangeToCheck] = determineRangeForGrid()
    
    % Initialization
    rangeToCheck = nan(nObjects,2);
    
    % Check provided sample points
    if obj.getShowData
        for iKrigingObj = 1:nObjects
            iObj = krigingObjVector(iKrigingObj);
            rangeToCheck(iKrigingObj,:) = [min(obj.KrigingObjects{iObj}.getOutputData)*-obj.MinMax(iObj),...
                                    max(obj.KrigingObjects{iObj}.getOutputData)*-obj.MinMax(iObj)];
        end
    end

    % Check conditional simulations
    for iRealizationNested = 1 : nRealizations
        for iKrigingObj = 1:nObjects
            if size(rangeToCheck,1)<iKrigingObj
                rangeToCheck(iKrigingObj,:) = [min(paretoSets{iRealizationNested}(:,iKrigingObj)),...
                                               max(paretoSets{iRealizationNested}(:,iKrigingObj))];
            else
                
                pointsToCheck = [rangeToCheck(iKrigingObj,1),...
                                 rangeToCheck(iKrigingObj,2),...
                                 (paretoSets{iRealizationNested}(:,iKrigingObj))'];
                rangeToCheck(iKrigingObj,:) = ...
                    [min(pointsToCheck),...
                     max(pointsToCheck)];
            end
        end
    end
    

    % Check if scale of data prediction are similar
    if obj.getShowData
        rangeToCheckProto = zeros(nObjects,2);
        for iKrigingObj = 1:nObjects
            rangeToCheckProto(iKrigingObj,:) = [min([rangeToCheckProto(iKrigingObj,1),rangeToCheckProto(iKrigingObj,2),min(paretoSets{iRealizationNested}(:,iKrigingObj))]),...
                                    max([rangeToCheckProto(iKrigingObj,1),rangeToCheckProto(iKrigingObj,2),min(paretoSets{iRealizationNested}(:,iKrigingObj))])];
        end
        deltaProto = rangeToCheckProto(:,2)-rangeToCheckProto(:,1);
        deltaCurrent = rangeToCheck(:,2)-rangeToCheck(:,1);
        if any((abs(deltaProto-deltaCurrent./deltaProto))>1e1)
            warning('Using data for interpolation range of solution space extends significantly the investigate soloution space. You might set obj.ShowData=false')
        end
    end
    
    if ~isempty(OutputConstraints)
        if ~all(size(rangeToCheck)==size(OutputConstraints))
            error('OutputConstraints should be of Size nKrigingObjectsX2')
        end
        rangeToCheckProto = bsxfun(@times,(-obj.MinMax(krigingObjVector)'),OutputConstraints);
        rangeToCheckProto = sort(rangeToCheckProto,2); % Guaranty that LBis first column
        
        rangeToCheck(~isinf(rangeToCheckProto)) = rangeToCheckProto(~isinf(rangeToCheckProto));
        
    end
end

% ----------------------------------------------------------
function [expectedParetoPoints] = determineExpectedParetoCurve()
    epsilon = 1e-10;
    lowerLevel = 0;
    upperLevel = 1;
    deltaNew = inf;
    deltaOld = inf;
    HV_QuantileOld = inf;
    upperLevelOld = upperLevel;
    lowerLevelOld = lowerLevel;
    
    firstTime = true;
    while upperLevel-lowerLevel>epsilon
        if ~firstTime
            deltaOld = deltaNew;
            HV_QuantileOld = HV_Quantile;
        end
        
        pointsSelected = pCover>=(upperLevel+lowerLevel)/2;
        HV_Quantile = Hypervolume_MEX(gridX(pointsSelected,:),bsxfun(@times,obj.ReferencePointHyperVolume(krigingObjVector),--obj.MinMax(krigingObjVector)));

        
        deltaNew = abs(HV_Quantile-HV_Mean);
        % Reverse process if the result is worse
        if (~firstTime)&&(deltaOld<deltaNew)
            upperLevel = upperLevelOld;
            lowerLevel = lowerLevelOld;
            if ~(HV_QuantileOld<HV_Mean)
                upperLevel = (upperLevelOld+lowerLevelOld)/2;
            else
                lowerLevel = (upperLevelOld+lowerLevelOld)/2;
            end
        else
            upperLevelOld = upperLevel;
            lowerLevelOld = lowerLevel;
            if HV_Quantile<HV_Mean
                upperLevel = (upperLevel+lowerLevel)/2;
            else
                lowerLevel = (upperLevel+lowerLevel)/2;
            end
        end
        
        firstTime=false;
    end
    
    HV_Quantile = Hypervolume_MEX(gridX(pointsSelected,:),bsxfun(@times,obj.ReferencePointHyperVolume(krigingObjVector),-obj.MinMax(krigingObjVector)));

    quantileNumber = (upperLevel+lowerLevel)/2;

    % Determine Expected Pareto Curve
    expectedDominatedPointsBool = pCover>=quantileNumber;

    if sum(expectedDominatedPointsBool)>0
        expectedParetoPoints = determineParetoSet_Mex(gridX(expectedDominatedPointsBool,:));
    else
        error('no point is dominated byexpected pareto curve')
    end
end

% ----------------------------------------------------------
function [sampleLocations] = createGrid()

    % Check if statistic toolbox can be used
    necessaryFeatures = {'Statistics_Toolbox'};
    validLincence = cellfun(@(f) license('checkout',f),necessaryFeatures);


    InputRange = [obj.LBInputVarInterpolation{krigingObjVector(1)}',...
                 obj.UBInputVarInterpolation{krigingObjVector(1)}'];

    if ~validLincence
        % Try to come as close as possible to the desired number of
        % samples
        nSamplePointsEachDim = round((nSampleLocations)^(1/nInputVar));
        nSamplePoints = nSamplePointsEachDim^(nInputVar);
        if abs(nSampleLocations-nSamplePoints)~=0
            warning('For creation of a grid (nSampleLocations)^(1/nInputVar) should be a natural number. nSampleLocations was corrected to %i',nSamplePoints)
        end

        % Assume that input space is the same for all objectives of
        % interest
        sampleLocationsProto = createNDGRID(InputRange(:,1),InputRange(:,2),nSamplePointsEachDim);

        % Randomize Locations
        sampleLocations = zeros(nSamplePoints,nInputVar*nRealizations);
        for iRealizationNested=1:nRealizations
            indicesSampleLocations=randi(nSamplePoints,nSamplePoints,1);
            sampleLocations(:,(iRealizationNested-1)*nInputVar+1:iRealizationNested*nInputVar) =...
                sampleLocationsProto(indicesSampleLocations,:);
        end
    else
        nSamplePoints = nSampleLocations;
        
        if obj.RepeatDesign
            nRealizationsNested = 1;
            % Is changed after loop again
            sampleLocations = zeros(nSamplePoints,nInputVar); 
        else
            nRealizationsNested = nRealizations;
            sampleLocations = zeros(nSamplePoints,nInputVar*nRealizations);
        end
        
        for iReal = 1: nRealizationsNested
            sampleLocationsProto = lhsdesign(nSamplePoints,nInputVar);
            sampleLocationsProto(:,:) = bsxfun(@times,sampleLocationsProto,InputRange(:,2)'-InputRange(:,1)');
            sampleLocationsProto(:,:) = bsxfun(@plus,sampleLocationsProto,InputRange(:,1)');
            sampleLocations(:,nInputVar*(iReal-1)+1:nInputVar*iReal) = sampleLocationsProto(:,:);
        end
        
        if obj.RepeatDesign
            sampleLocations = repmat(sampleLocationsProto,1,nRealizations);
        end

    end

    if ~isempty(inequalityConstraintInput)
        sampleLocations = sampleLocations(inequalityConstraintInput(sampleLocations),:);
        nSampleLocations = size(sampleLocations,1);
    end

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
