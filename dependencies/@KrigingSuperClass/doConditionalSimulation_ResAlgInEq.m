function [realizationCurve,ReturnSamplePoints] = doConditionalSimulation_ResAlgInEq(obj,varargin)
% Either:
% [realizationCurve,ReturnSamplePoints] = doConditionalSimulation_ResAlgInEq(sampleLocations,inEqConstraint)
% Or: 
% [realizationCurve,ReturnSamplePoints] = doConditionalSimulation_ResAlgInEq(InputRange,nRealization,accuracy,inEqConstraint)
%
% Conditional Simulation are tools for evaluate the impact of uncertainity
% on the results of compex procedures. They are so-called Monte Carlo
% methods. For further information see (Chiles 1999 - "Geostatistics:
% modeling spatial uncertainty") 
%
% This function uses "Kriging residual algorithm" (see chevalier2015 - Fast
% update of conditional simulation ensembles):
% 1. NonConditional Simulation: Generating realizations of a process with
%    kriging mean and in agreement with the Kriging Covarisance matrix
%    z_noncond(x)~N(m(x),C(x)) 
% 2. Conditioning on the Data: z_cond(x) = m(x) + r(x)
%       m(x) ... kriging estimation
%       r(x) ... residual estimation: r(x) = z_noncond(x) -
%                                            w(x)'z_noncond(x_input)
%       w(x) ... kriging weights used for calculating m(x)
%       z_noncond(x_input) ... nonconditional simulation at the user  
%                              provided input data points.
%    I.e., the residual represents the modeling error by replacing the
%    actual output data with the nonconditional simulation values at the
%    particular points
%
% Input: 
% sampleLocations ... contains the simulation points at which conditional
%                     realization shall be drawn from. The finer the grid
%                     the smoother the realization curve.
%                     (nSamplePointsXnInputVar*nRealizations)
% InputRange ... contains the min/max values in the first and second
%                column, respectively (nInputVarX2)
% nRealization ... number of realization curves which shall be drawn
%                  (scalar)
% accuracy ... number of points in each direction in the grid. Based on
%              that the "marginal" width between points is ~1/accuracy
%              (Scalar)
% inEqConstraint ... handle to a function of the form 
%                    "function [conditionHoldsBool] = inEqConstraint(simulationOutput)"
%                    Gives out a boolean vector containing true at the
%                    entries where the rows of the matrix simulationOutput
%                    do not hold the inequality constraints
%
% Output: 
% realizationCurve ... simulated realization curves
%                      (nSamplePointsXnRealization)
% ReturnSamplePoints ... simualtion points at which the realizations are
%                        drawn (nSamplePointsXnInputVar)
%
% You can set: 
% - ShowWaitingBar ... if true, the process will be displayed i na waiting
%                      bar
%
% You can get: -
%
% This function runs considerable faster when all realization are drawn at
% identical simulation points given by sampleLocations;
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

if obj.getShowWaitingBar
    hWaitBar = waitbar(0,'Conditional Simulation is Running');
    tic
end
%% Initialization

% Define sample locations for simulation,number of simulation points, and
% number of realizations
if length(varargin)==1||length(varargin)==2
    [sampleLocations,nSamplePoints,nRealization] = checkGivenMatrix;
else
    [sampleLocations,nSamplePoints,nRealization] = checkGivenInputRange;
end

% Normalize Input if necessary
if obj.NormInput==1
    % Save original data for later
    nonNormSampleLocations = sampleLocations;
    
    % Scale sample points of each realization to the range [0,1]
    for iRealization=1:nRealization
        for iInput = 1:obj.nInputVar
            column = iInput +(iRealization-1)*obj.nInputVar;
            sampleLocations(:,column) = obj.scale(sampleLocations(:,column), min(obj.InputData_True(:,iInput)),max(obj.InputData_True(:,iInput)),0,1);
        end
    end
else
    nonNormSampleLocations = sampleLocations;
end 


% Backup in order to return to old state after conditional simulation
maxSizeOfPredictionsBackup = obj.getMaxSizeOfPredictions;
UseGPRMatlabBackup = obj.getUseMatlabRegressionGP;
% HeterogeneousNoiseBackUp = obj.getHeterogeneousNoise;
obj.setUseMatlabRegressionGP(false)

% Allocate Memory: 
    % Save the resulting realization curves
realizationCurve = zeros(nSamplePoints,nRealization);
    % Used for saving the estimatied model residual in step 2
randomResidual = zeros(nSamplePoints,1);

% Make sure that no permanent parameter changes happens when porogram
% crashes
% try
%% Step 0: Generate normal Kriging estimation at simulation points and save associated kriging weights

    % Do predictions in one step
    obj.setMaxSizeOfPredictions(obj.getnExperiments+nSamplePoints)

    % Actual Prediction. Do not allow heterogenous noise here as it is not
    % used for non-conditional simulation
    outputMatrix = obj.prediction(nonNormSampleLocations(:,1:obj.getnInputVar));
    
    % Save results
    krigingWeights = obj.getWeights;
    krigingWeights = krigingWeights(1:obj.getnExperiments,:);
    meanEstimation = outputMatrix(:,1);
    if obj.NormOutput
        meanEstimation = obj.scale(meanEstimation,...
                                   min(obj.getOutputData),max(obj.getOutputData),...
                                   0,1);
    end
%     stdEstimation = outputMatrix(:,2); 
    
    % Intiital evaluation 
    processOld = 0;

%% Step 1: Do nonconditional simulation at data points

    % Check if simulation points of all relization are identical. If yes, non
    % conditional simulation can be done once in a parallel fashion
    allRealizationPointsAreTheSame = true;
    for iRealization = 2:nRealization
        columnSamplePoints1 = 1:1*obj.nInputVar;
        columnSamplePoints2 = (iRealization-1)*obj.nInputVar+1:iRealization*obj.nInputVar;
        allRealizationPointsAreTheSame = allRealizationPointsAreTheSame&...
            sum(sum(abs(sampleLocations(:,columnSamplePoints1)-sampleLocations(:,columnSamplePoints2))))==0;
    end

    % Do Actual nonconditional simulation
    if allRealizationPointsAreTheSame
        
        % Calculated Kriging covariance matrix w.r.t to all simulation
        % points
        columnSamplePoints = 1:1*obj.nInputVar;
        covMatrix = createCovMatrix(sampleLocations(:,columnSamplePoints));
        
        indices1 = 1:size(covMatrix,1)+1:size(covMatrix,1)*obj.nExperiments;
        indices2 = indices1(end) + size(covMatrix,1)+1:size(covMatrix,1)+1:size(covMatrix,1)^2;

        if (isempty(obj.KriKitObjNoise)||isempty(obj.KriKitObjNoise.getOutputData))
            covMatrix([indices1,indices2]) = covMatrix([indices1,indices2]) + (obj.sigmaError).^2;    
        else
            predLogNoise = obj.KriKitObjNoise.prediction(nonNormSampleLocations(:,1:obj.getnInputVar));
            covMatrix([indices1,indices2]) = bsxfun(@plus,covMatrix([indices1,indices2]), exp(predLogNoise(:,1))' );
        end
        
        % Random drawing following the distribution ~N(0,covMatrix)
        % First nDataXnData part of the covMatrix has to contain the data
        % points covariance matrix, where sigmaError is added to the
        % diagonal entries (Noisy measurements)
        [unconditionedReal] = generatedNonconditionalSImulation(covMatrix,nRealization);
    end

%% Step 2: Conditional simulation
    if allRealizationPointsAreTheSame
        randomResidual = zeros(nSamplePoints,nRealization);
        randomResidual(:) = unconditionedReal(:,obj.getnExperiments+1:end)'- krigingWeights(:,obj.getnExperiments+1:end)'*unconditionedReal(:,1:obj.getnExperiments)';
        realizationCurve = bsxfun(@plus,randomResidual,meanEstimation(obj.getnExperiments+1:end));
    else

        % Has to be done individual for each realization
        for iRealization = 1:nRealization
            % Do nonconditional simulation (Only needed when simulation point
            % differ between realizations
            columnSamplePoints = (iRealization-1)*obj.nInputVar+1:iRealization*obj.nInputVar;
            
            % Actual Prediction
            outputMatrix = obj.prediction(nonNormSampleLocations(:,columnSamplePoints));


            % Save results
            krigingWeights = obj.getWeights;
            krigingWeights = krigingWeights(1:obj.getnExperiments,:);
            meanEstimation = outputMatrix(:,1);
            covMatrix = createCovMatrix(sampleLocations(:,columnSamplePoints));
            if ~isempty(obj.HeterogeneousNoise)
                covMatrix(1:size(covMatrix,1)+1:size(covMatrix,1)*obj.nExperiments) = covMatrix(1:size(covMatrix,1)+1:size(covMatrix,1)*obj.nExperiments) + (obj.HeterogeneousNoise').^2;
            else
                covMatrix(1:size(covMatrix,1)+1:size(covMatrix,1)*obj.nExperiments) = covMatrix(1:size(covMatrix,1)+1:size(covMatrix,1)*obj.nExperiments) + (obj.sigmaError).^2;
            end
             
            [unconditionedReal] = generatedNonconditionalSImulation(covMatrix,nRealization);

            % Do Conditional simulation constrainted by inEqConstraint
            realizationFound = false;
            while ~realizationFound

                % Calculate estimation modeling error
                randomResidual(:) = unconditionedReal(iRealization,obj.getnExperiments+1:end)'- krigingWeights(:,obj.getnExperiments+1:end)'*unconditionedReal(iRealization,1:obj.getnExperiments)';

                % Final calculation of conditional simulation
                realizationCurve(:,iRealization) = meanEstimation(obj.getnExperiments+1:end) + randomResidual;

                % Update Process bar if wanted
                if obj.getShowWaitingBar
                    processNew = iRealization/nRealization;
                    if (processNew-processOld)>0.01
                        waitbar(processNew,hWaitBar,sprintf('Process For Conditional Simulation: %f',processNew))
                        processOld = processNew;
                    end
                end
            end


        end
    end
%% Final things
    if obj.NormOutput
        realizationCurve = obj.scale(realizationCurve,0,1,...
                                   min(obj.getOutputData),max(obj.getOutputData));
    end
    
    % Define output
    ReturnSamplePoints = nonNormSampleLocations(obj.nExperiments+1:end,:);

    % Close the waiting bar
    if obj.getShowWaitingBar
        close(hWaitBar)
    end

    obj.setUseMatlabRegressionGP(UseGPRMatlabBackup);
    obj.setMaxSizeOfPredictions(maxSizeOfPredictionsBackup);

%% Nested Functions 
% -------------------------------------------------------------------------
function [unconditionedReal] = generatedNonconditionalSImulation(covMatrix,nRealization)
    necessaryFeatures = {'Statistics_Toolbox'};
    validLincence = cellfun(@(f) license('checkout',f),necessaryFeatures);
    if validLincence
        try
            unconditionedReal= mvnrnd(meanEstimation,covMatrix,nRealization);
        catch ex
            warning(ex.message);
            covMatrix = correctCovMatrix(covMatrix);
            unconditionedReal= mvnrnd(meanEstimation,covMatrix,nRealization);
        end
    else
        [eigVector,eigValue]=eig(covMatrix);
        if any(eigValue<0)
            covMatrix = correctCovMatrix(covMatrix);
            [eigVector,eigValue]=eig(covMatrix);
        end
        unconditionedReal = randn(nRealization,obj.getnExperiments+nSamplePoints);
        unconditionedReal = eigVector*sqrt(eigValue)*eigVector'*unconditionedReal';
        unconditionedReal =  bsxfun(@plus,unconditionedReal,meanEstimation);
        unconditionedReal = unconditionedReal';
    end
    
end
% -------------------------------------------------------------------------
function [covMatrix] = correctCovMatrix(covMatrix)
    % Covariance matrix has to be semipositive definite
    eigBefore = eig(covMatrix);
    nRowsColumns = size(covMatrix,1);
    indicesDiag = 1:nRowsColumns+1:nRowsColumns^2;
    covMatrix(indicesDiag) = covMatrix(indicesDiag)+covMatrix(indicesDiag)*1e-10;
%             covMatrix(indicesDiag) = covMatrix(indicesDiag)+min(eigBefore);
    eigAfter = eig(covMatrix);
    warning('covMatrix is not positive semidefinite (smallest eigenvalue is %g). \nAfter adding 1e-10 to diagonal, max Diff in eigenvalues is : %g\n',min(eigBefore),max(abs(eigBefore-eigAfter)))
end
% -------------------------------------------------------------------------
function [covMatrix] = createCovMatrix(sampleLocations)
    nPoints = nSamplePoints+obj.nExperiments;
    
    [a,b] = ndgrid(1:nPoints,1:nPoints);
    combinationMatrix = [b(:),a(:)];
    if(obj.CovariogramUsesEuclideanDistance)
        delta = sqrt(sum((bsxfun(@minus,sampleLocations(combinationMatrix(:,1),:),sampleLocations(combinationMatrix(:,2),:)).^2),2));
        nRowsColumns = sqrt(size(delta,1));
        delta = reshape(delta,nRowsColumns,nRowsColumns);

        % Calculate the variogram values for theses distances
        covMatrix =  obj.CovarModel(delta,0);
        
    elseif(obj.CovariogramUsesAbsoluteDistance)
        
        % Calculate the covariogram values for absolute input distances 
        delta=abs(bsxfun(@minus,sampleLocations(combinationMatrix(:,1),:),sampleLocations(combinationMatrix(:,2),:)));
        nRowsColumns=sqrt(size(delta,1));
        d = ones(nRowsColumns,nRowsColumns);
        for iN = 1:nRowsColumns
            for iV = 1:obj.nInputVar
               d(:,(iN-1)*obj.nInputVar+iV)=delta(nRowsColumns*(iN-1)+1:nRowsColumns*(iN),iV);
            end
        end
        delta = reshape(d,nRowsColumns,obj.nInputVar,nRowsColumns);
        covMatrix =  reshape(obj.CovarModel(delta,0),nRowsColumns,nRowsColumns);
        
    else
        error('CovariogramUsesEuclideanDistance or CovariogramUsesAbsoluteDistance should be defined')
    end
    
end
% -------------------------------------------------------------------------
function [sampleLocations,nSamplePoints,nRealization] = checkGivenMatrix
    sampleLocations = varargin{1};
    if mod(size(sampleLocations,2),obj.nInputVar)~=0||isempty(sampleLocations)
        error('"sampleLocations" has to be size of nSamplePointsX(nRealization*nInputVar(=%i))',obj.nInputVar)
    end

    nRealization = size(sampleLocations,2)/obj.nInputVar;

    nSamplePoints = size(sampleLocations,1);

    % Add Known Sample Locations where measurement are already made
    sampleLocations = [repmat(obj.getInputData,1,nRealization);sampleLocations];
end
% -------------------------------------------------------------------------
function [sampleLocations,nSamplePoints,nRealization] = checkGivenInputRange
    % Check Input
    InputRange = varargin{1};

    if length(varargin)>1&&~isempty(varargin{2})
        nRealization = varargin{2};
    else
        nRealization = 1;
    end

    if length(varargin)>2&&~isempty(varargin{3})
        accuracy = varargin{3};
    else
        accuracy = 100;
    end

    % Create Sample Points
    nSamplePoints = accuracy^obj.nInputVar;
    necessaryFeatures = {'Statistics_Toolbox'};
    validLincence = cellfun(@(f) license('checkout',f),necessaryFeatures);

    if ~validLincence
        sampleLocationsProto = createNDGRID(InputRange(:,1),InputRange(:,2),accuracy);

        % Randomize Locations
        sampleLocations = zeros(nSamplePoints+obj.nExperiments,obj.nInputVar*nRealization);
        for iRealizationNested=1:nRealization
            indicesSampleLocations=randi(nSamplePoints,nSamplePoints,1);
            sampleLocations(:,(iRealizationNested-1)*obj.nInputVar+1:iRealizationNested*obj.nInputVar) = [obj.getInputData;sampleLocationsProto(indicesSampleLocations,:)];
        end
    else
        sampleLocations = zeros(nSamplePoints+obj.nExperiments,obj.nInputVar*nRealization);
        for iRealizationNested=1:nRealization
            sampleLocationsProto = lhsdesign(nSamplePoints,obj.nInputVar);
            sampleLocationsProto(:,:) = bsxfun(@times,sampleLocationsProto,InputRange(:,2)'-InputRange(:,1)');
            sampleLocations(:,(iRealizationNested-1)*obj.nInputVar+1:iRealizationNested*obj.nInputVar) =...
                [obj.getInputData;bsxfun(@plus,sampleLocationsProto,InputRange(:,1)')];
        end
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
