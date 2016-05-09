function [realizationCurve,ReturnSamplePoints] = doConditionalSimulation_SeqSim(obj,varargin)
% Either:
% [realizationCurve,ReturnSamplePoints] = doConditionalSimulation_SeqSim(sampleLocations)
% Or: 
% [realizationCurve,ReturnSamplePoints] = doConditionalSimulation_SeqSim(InputRange,nRealization,accuracy)
%
% Conditional Simulation are tools for evaluate the impact of uncertainity
% on the results of compex procedures. They are so-called Monte Carlo
% methods. For further information see (Chiles 1999 - "Geostatistics:
% modeling spatial uncertainty") 
%
% This function uses "Sequential Simulation". I.e., realizations at each
% simulation points are sequentially calculated. Once a realization
% calculated, it is added to the current data set. The kriging model is
% then updated (including covariogram parameter estimation and calculation
% of the covariogram matrix and its inverse). NOTE: Sequential simulation
% is considerably slower than the"Kriging residual algorithm" (see
% doConditionalSimulation_ResAlgInEq)
%
% Input: 
% sampleLocations ... contains the simulation points at which conditional
%                     realization shall be drawn from. THe finer the grid
%                     the smoother the realization curve.
%                     (nSamplePointsXnInputVar)
% InputRange ... contains the min/max values in the first and second
%                column, respectively (nInputVarX2)
% nRealization ... number of realization curves which shall be drawn
%                  (scalar)
% accuracy ... number of points in each direction in the grid. Based on
%              that the "marginal" width between points is ~1/accuracy
%              (Scalar)
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
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

useSparse = false; 

if obj.getShowWaitingBar
    hWaitBar = waitbar(0,'Conditional Simulation is Running');
    tic
end
%% initialization
	if length(varargin)==1
        [sampleLocations,nSamplePoints,nRealization] = checkGivenMatrix;
        
    else
        [sampleLocations,nSamplePoints,nRealization] = checkGivenInputRange;
    end
    
    % Normalize Input if necessary
    if obj.NormInput==1
        nonNormSampleLocations = sampleLocations;
        
        for iRealization=1:nRealization
            for iInput = 1:obj.nInputVar
                column = iInput +(iRealization-1)*obj.nInputVar;
                sampleLocations(:,column) = obj.scale(sampleLocations(:,column), min(obj.InputData_True(:,iInput)),max(obj.InputData_True(:,iInput)),0,1);
            end
        end
    else
        nonNormSampleLocations = sampleLocations;
    end
%% 
        
        % Allocate Memory
            % Covariance Matrix for the measurements
        obj.calcCovariogramMatrix
        CovData = obj.CovariogramMatrix(1:obj.nExperiments,1:obj.nExperiments);
            % Evaluation of the basis functions at the measurement points
        basisFunctionEval = obj.CovariogramMatrix(obj.nExperiments+1:end,1:obj.nExperiments);
            % Matrix for saving distance between one sample point and the
            % measurements. Used for the calculation of the covariance
        distanceVec = zeros(obj.nExperiments+nSamplePoints,1);
            % Matrix for saving covariance between one sample point and the
            % measurements. 
        newCovVec = zeros(obj.nExperiments+nSamplePoints,1);
            % Save sample points which are added to the sample space
            % succesively 
        inputValues = zeros(obj.nExperiments+nSamplePoints,obj.nInputVar);
        outputValues = zeros(obj.nExperiments+nSamplePoints,1);
            % Save the resulting realization curves
        realizationCurve = zeros(nSamplePoints,nRealization);
            % Random number following a normal distribution
        randomVec = randn(nSamplePoints,nRealization);
        % 
        
        
        processOld = 0;
        for iRealization = 1:nRealization
            columnSamplePoints = (iRealization-1)*obj.nInputVar+1:iRealization*obj.nInputVar;
            
            % Initialize Covariance matrix 
            if useSparse
                CovDatatotal = sparse(zeros(obj.nExperiments+nSamplePoints+obj.nBasisFct,...
                             obj.nExperiments+nSamplePoints+obj.nBasisFct));
            else
                CovDatatotal = zeros(obj.nExperiments+nSamplePoints+obj.nBasisFct,...
                             obj.nExperiments+nSamplePoints+obj.nBasisFct);
                invCovDatatotal = zeros(obj.nExperiments+nSamplePoints+obj.nBasisFct,...
                             obj.nExperiments+nSamplePoints+obj.nBasisFct);
                invCovDataOnly = zeros(obj.nExperiments+nSamplePoints,...
                             obj.nExperiments+nSamplePoints);
                basisMatrixTotal = zeros(obj.nBasisFct,obj.nExperiments+nSamplePoints);
            end
            CovDatatotal(1:obj.nExperiments,1:obj.nExperiments) = CovData;
            CovDatatotal(obj.nExperiments+1:obj.nExperiments+obj.nBasisFct,1:obj.nExperiments) = basisFunctionEval;
            CovDatatotal(1:obj.nExperiments,obj.nExperiments+1:obj.nExperiments+obj.nBasisFct) = basisFunctionEval';
            basisMatrixTotal(:,1:obj.nExperiments) = basisFunctionEval;
            
            for iSamplePoint=1:nSamplePoints
                
                if obj.getShowWaitingBar
                    processNew = ((iRealization-1)*nSamplePoints + iSamplePoint)/(nRealization*nSamplePoints);
                    if (processNew-processOld)>0.01
                        waitbar(processNew,hWaitBar,sprintf('Process For Conditional Simulation: %f',processNew))
%                         toc
                        processOld = processNew;
                    end
                    
                end
                
                indexTestValue = obj.nExperiments+iSamplePoint-1;
                
                % Update Sample Points which are already added to sample
                % space
                inputValues(1:indexTestValue,:) = sampleLocations(1:indexTestValue,columnSamplePoints);
                outputValues(1:indexTestValue,:) = [obj.OutputData;realizationCurve(1:iSamplePoint-1,iRealization)];

                % calculate covariance between new sample point and already
                % added sample point
                if obj.CovariogramUsesEuclideanDistance
                    newCovVec(1:indexTestValue)=calcCovargramVectorEuclideanDistance(inputValues(1:indexTestValue,:),sampleLocations(indexTestValue+1,columnSamplePoints));
                elseif obj.CovariogramUsesAbsoluteDistance
                    newCovVec(1:indexTestValue)=calcCovargramVectorAbsoluteDistance(inputValues(1:indexTestValue,:),sampleLocations(indexTestValue+1,columnSamplePoints));
                else
                    error('EitherCovariogramUsesEuclideanDistance or CovariogramUsesAbsoluteDistance have to be set true')
                end
                
                % Final Covariogram Vecotr
                % Extend by Basis Functions
                for iBasis1 = 1 : obj.nBasisFct
                    basis = obj.BasisFct{iBasis1}(obj.BasisFctParameters,nonNormSampleLocations(indexTestValue+1,columnSamplePoints));
                    newCovVec(indexTestValue+iBasis1) = basis;
                end
                
                % Calculate Actual Prediction
                if obj.nBasisFctParameters == 0&&obj.getUseSimpleKriging==1
                    if iSamplePoint==1
                        invCovDatatotal(1:indexTestValue,1:indexTestValue) =...
                            inv(CovDatatotal(1:indexTestValue,1:indexTestValue));
                        invCovDataOnly(1:indexTestValue,1:indexTestValue) = invCovDatatotal(1:indexTestValue,1:indexTestValue);
                    end
%                     krigingWeights = CovDatatotal(1:indexTestValue,1:indexTestValue)\newCovVec(1:indexTestValue,:);
                    krigingWeights = invCovDatatotal(1:indexTestValue,1:indexTestValue)*newCovVec(1:indexTestValue,:);
                else
                    if iSamplePoint==1
                        invCovDatatotal(1:indexTestValue+obj.nBasisFct,1:indexTestValue+obj.nBasisFct) =...
                            inv(CovDatatotal(1:indexTestValue+obj.nBasisFct,1:indexTestValue+obj.nBasisFct));
                        invCovDataOnly(1:indexTestValue,1:indexTestValue) = inv(CovDatatotal(1:indexTestValue,1:indexTestValue));
                    end
                    
                    % Use gaussian elemination instead of the inverting the matrix
%                     krigingWeightsCheck = CovDatatotal(1:indexTestValue+obj.nBasisFct,1:indexTestValue+obj.nBasisFct)\newCovVec(1:indexTestValue+obj.nBasisFct,:);
                    krigingWeights = invCovDatatotal(1:indexTestValue+obj.nBasisFct,1:indexTestValue+obj.nBasisFct)*newCovVec(1:indexTestValue+obj.nBasisFct,:);
                    
%                     sum(abs(krigingWeightsCheck-krigingWeights))
                end
                
                
                
                % Final Prediction
                prediction = [outputValues(1:indexTestValue,:);zeros(size(krigingWeights,1)-indexTestValue,1)]'*krigingWeights;
                sigmaEstimation=calcSigmaEstimation(krigingWeights);
                
                
                if obj.NormOutput
                    % Convert just the mean to original scale (not sigmaEstimation!)
                    prediction = obj.scale(prediction,0,1,...
                                       min(obj.getOutputData),max(obj.getOutputData));
                    
                    realizationCurve(iSamplePoint,iRealization) = randomVec(iSamplePoint,iRealization)*sigmaEstimation + prediction;
                    
                    % Convert result to normalized scale (needed for remaining iterations)
                    realizationCurve(iSamplePoint,iRealization) = obj.scale(realizationCurve(iSamplePoint,iRealization),...
                                       min(obj.getOutputData),max(obj.getOutputData),0,1);
                       
                else
                    realizationCurve(iSamplePoint,iRealization) = randomVec(iSamplePoint,iRealization)*sigmaEstimation + prediction;
                end
                
                
                % Update CovDatatotal
                    % For better overview defined indices
                indicesOld = 1:indexTestValue;
                indicesNew = indexTestValue+1;
                indicesNewBasis = indexTestValue+2:indexTestValue+1+obj.nBasisFct;
                indicesOldBasis = indexTestValue+1:indexTestValue+obj.nBasisFct;
                
                
                % !!!! I do not know why this part is important, but
                % strange things are happening as soon as this part is
                % not existing anymore ..
                    % Update Basis function entries (first shift old
                    % entries
                CovDatatotal(indicesNewBasis,indicesOld) = CovDatatotal(indicesOldBasis,indicesOld);
                CovDatatotal(indicesOld,indicesNewBasis) = CovDatatotal(indicesOld,indicesOldBasis);
                    % Add new entry for Basis fucntion evaluation
                CovDatatotal(indicesNewBasis,indicesNew) = newCovVec(indicesOldBasis);
                CovDatatotal(indicesNew,indicesNewBasis) = newCovVec(indicesOldBasis)';
                    % Update Covariance between samples points 
                CovDatatotal(indicesNew,indicesOld) = newCovVec(indicesOld)';
                CovDatatotal(indicesOld,indicesNew) = newCovVec(indicesOld);
                    % Update Diagonal (just extending)
                CovDatatotal(indicesNew,indicesNew) = CovDatatotal(1,1);

                
                invCovDataOnly(1:indexTestValue+1,1:indexTestValue+1) = ...
                    invupdateapp(invCovDataOnly(1:indexTestValue,1:indexTestValue),newCovVec(1:indexTestValue),newCovVec(1:indexTestValue)',CovDatatotal(1,1));
                invCovDatatotal(1:indexTestValue+1,1:indexTestValue+1) = invCovDataOnly(1:indexTestValue+1,1:indexTestValue+1);
                
                if obj.nBasisFct>=1
                    basisMatrixTotal(:,indexTestValue+1) = newCovVec(indexTestValue+1:indexTestValue+obj.nBasisFct);
                end
                for iBasis1=1:obj.nBasisFct
                    addVec = [basisMatrixTotal(iBasis1,1:indexTestValue+1),zeros(1,iBasis1-1)];
                    invCovDatatotal(1:indexTestValue+1+iBasis1,1:indexTestValue+1+iBasis1) = ...
                        invupdateapp(invCovDatatotal(1:indexTestValue+1+iBasis1-1,1:indexTestValue+1+iBasis1-1),addVec',addVec,0);
                end
                
               
            end
        end
        %%

        if obj.NormOutput
            for iRealization = 1:nRealization
                realizationCurve(:,iRealization) = obj.scale(realizationCurve(:,iRealization),0,1,...
                               min(obj.getOutputData),max(obj.getOutputData));
                % Do not normalize standard deviation as for the calculation no
                % normalized output values are used!!!!!!!!!!!!!!!!!!
            end
        end
        ReturnSamplePoints = nonNormSampleLocations(obj.nExperiments+1:end,:);

        % Close the waiting bar
        if obj.getShowWaitingBar
            close(hWaitBar)
        end
        
%% Nested Functions 
    function [CovargramVector] = calcCovargramVectorEuclideanDistance(upToNowSamplePoints,newSamplePoint)
            % Calculate actual difference
        distance = bsxfun(@minus,upToNowSamplePoints,newSamplePoint);
            % Euclidean norm
        distance = sqrt(sum(distance.^2,2));
        % Calculate variogram vector
            % Use calculated distances for the determination of the
            % variogram vector
        CovargramVector(1:indexTestValue,:) = obj.CovarModel(distance,obj.AbsoluteInterpolation);
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
% -------------------------------------------------------------------------
    function [CovargramVector] = calcCovargramVectorAbsoluteDistance(upToNowSamplePoints,newSamplePoint)
            % Calculate actual difference
            distance = bsxfun(@minus,upToNowSamplePoints,newSamplePoint);
                % Absolute Value
            distance = abs(distance);
            % Calculate variogram vector
                % Use calculated distances for the determination of the
                % variogram vector
            CovargramVector(1:indexTestValue,:) = obj.CovarModel(distance,obj.AbsoluteInterpolation);
    end
% -------------------------------------------------------------------------
function [sigmaEstimation]=calcSigmaEstimation(krigingWeights)
% For a better commented version see "prediction"

        covariance_zero  = obj.CovarModel(zeros(1,size(obj.DistInput,2)),1);
        sigmaEstimation = covariance_zero - 2*krigingWeights'*newCovVec(1:indexTestValue+obj.nBasisFct,:) +...
                          krigingWeights'* CovDatatotal(1:indexTestValue+obj.nBasisFct,1:indexTestValue+obj.nBasisFct)*krigingWeights;

        % Only the diagonal are the number we are looking for
        sigmaEstimation = diag(sigmaEstimation);
        sigmaEstimation = (sigmaEstimation).^(1/2);
        
        % Check for imaginary numbers
        testImag = find(imag(sigmaEstimation)~=0, 1);
        if(~isempty(testImag))
            sigmaEstimation = real(sigmaEstimation);
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
