function [OutputTotal] = prediction(obj,input)
% [Output] = prediction(input)
% This function predicts the Outpus Values for the Given Input Points
% input ... data points where the function value shall be predicted (nxInputVar Matrix)
% Output: Is a matrix nX2. n is the number of given data points
%         First column contains the actual prediction 
%         Second column contains the estimated prediction error
%
% You can set:
% - AbsoluteInterpolation ... Force estimation to pass all sample points
% - MaxSizeOfPredictions ... % Prediction gets very slow when the input
%                              matrix is getting too big.
%                              MaxSizeOfPredictions is the maximal number
%                              of points at which the prediction is
%                              parallel calculated. By default 100
% - ShowWaitingBar ... if true, process is displayed in a waitbar 
% - UseMatlabRegressionGP ... if true MathWorks Gaussian Process Regression
%                             function is used (Matlab2015b or newer and
%                             statistic toolbox are needed)
% - BasisFct ... define basis function for better Kriging prediction
% - UseSimpleKriging ... Simple Kriging: Linear coefficients of the basis
%                        function is not automaically estimated and is set
%                        to 1. 
%                        So only constant function can be considered
%                        "@(p,x)a". Default: false
%
% You can get:
% - Weights ... kriging wights used for the caluclation of "Output"
%               - krigingWeights = obj.InvCovariogramMatrix*obj.CovargramVectors;
%               - prediction = [obj.OutputData;zeros(nBasisFct,1)]'*krigingWeights;
%               Note: last nBasisFct entries, represent the laplace
%               multipliers for the basis function
%               
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.
    

    
    
    %% Check Conditions
    % Check if the Equation System is not singular
    if (obj.nExperiments<=obj.nBasisFct)
        error('Not enough sample Points to get a unique solution for the krigingWeights. \n The number of experiments/sample poins(here %i) should be bigger than the number of basis functions (here %i)',obj.nExperiments,obj.nBasisFct);
    end
    if size(input,2)~=obj.nInputVar&&size(input,1)==obj.nInputVar
        input = input';
        warning('Input vector should be a nx%i vector but it is a %ix%i vector. The vector is transposed',obj.nInputVar,size(input,1),size(input,2))
    end
    if size(input,2)~=obj.nInputVar&&size(input,1)~=obj.nInputVar
        error('Input vector should be a nx%i vector but it is a %ix%i vector',obj.nInputVar,size(input,1),size(input,2))
    end
    
    %% Prediction GPR
    % If MathWorks Gaussian Process Regression function is used
    if obj.UseMatlabRegressionGP&&(isempty(obj.KriKitObjNoise)||isempty(obj.KriKitObjNoise.getOutputData))
        if obj.NormInput==1
            for iInput = 1:obj.nInputVar
                input(:,iInput) = obj.scale(input(:,iInput), min(obj.InputData_True(:,iInput)),max(obj.InputData_True(:,iInput)),0,1);
            end
        end
        [prediction,sigmaEstimation] = predict(obj.GPR_Model,input);
        OutputTotal = [prediction,sigmaEstimation];
        
        if obj.NormOutput
            % Scale = (x-minOld)*((maxNew-minNew)/(maxOld-minOld)) + minNew
            %       = (x-a)*b + c
            % -> x~N((y-a)*b+c,sigma^2*b^2)
            OutputTotal(:,1) = obj.scale(OutputTotal(:,1),0,1,...
                           min(obj.getOutputData),max(obj.getOutputData));
            b = (max(obj.getOutputData)-min(obj.getOutputData))/(1-0);
            OutputTotal(:,2) = OutputTotal(:,2)*b;
            % Do not scale standard deviation as for the calculation no
            % normalized output values are used!!!!!!!!!!!!!!!!!!
            % f=(x-minX)./(maxX-minX).*(UB-LB)+LB;
        end
        
        return
    end
    
    %% Prediction KriKit
    % You may need the non-normalized input 
    nonNormInput = input;
    
    % Separate the prediction to avoid memory overload
    if size(input,1)>obj.MaxSizeOfPredictions
        
        % Initialization
        TotalInput = input;
        TotalNonNormInput =nonNormInput;
        OutputTotal = zeros(size(TotalInput,1),2);
%         obj.Weights = zeros(obj.getnExperiments+obj.getnBasisFct,size(TotalInput,1));
%         CovVecTot = zeros(obj.getnExperiments+obj.getnBasisFct,size(TotalInput,1));
        
        if obj.ShowWaitingBar==1;
            str =  horzcat('Please wait prediction runs ... already done: ',num2str(0),'%');
            hWaitingBar = waitbar(0,str);
        end
        predictionSteps = floor(size(TotalInput,1)/obj.MaxSizeOfPredictions);
            
        % Calculate the prediction for each part successively
        for iPredictions = 1 : predictionSteps;
            
            % Define current input
            indexBeginn = (iPredictions-1)*obj.MaxSizeOfPredictions+1;
            indexEnd = (iPredictions)*obj.MaxSizeOfPredictions;
            input = TotalInput(indexBeginn:indexEnd,:);
            nonNormInput = TotalNonNormInput(indexBeginn:indexEnd,:);
            nInput = size(input,1);
            
            % Nomalize Input (Cannot be separated from this place otherwise memeroy overload is risked)
            if obj.NormInput==1
                for iInput = 1:obj.nInputVar
                    input(:,iInput) = obj.scale(input(:,iInput), min(obj.InputData_True(:,iInput)),max(obj.InputData_True(:,iInput)),0,1);
                end
            end
            
            % Calculate the associated prediction (Output is calculated)
            doPredictionUniversalKriging();
            
            OutputTotal(indexBeginn:indexEnd,:) = Output;
            
            % Give Feedback
            if obj.ShowWaitingBar==1;
                str =  horzcat('Please wait prediction runs ... already done: ',num2str(iPredictions/predictionSteps*100,'%3.2f'),'%');
                waitbar(iPredictions/predictionSteps,hWaitingBar,str)
            end
        end
        
        % Remaining prediction points (remaining by devision with
        % MaxSizeOfPredictions) 
        input = TotalInput(indexEnd+1:end,:);
        nonNormInput = TotalNonNormInput(indexEnd+1:end,:);
        if obj.NormInput==1
            for iInput = 1:obj.nInputVar
                input(:,iInput) = obj.scale(input(:,iInput), min(obj.InputData_True(:,iInput)),max(obj.InputData_True(:,iInput)),0,1);
            end
        end
        nInput = size(input,1);
        doPredictionUniversalKriging;
        OutputTotal(indexEnd+1:end,:) = Output;
%         obj.Weights(:,indexEnd+1:end) = krigingWeights;
%         CovVecTot(:,indexEnd+1:end) = obj.CovargramVectors;
        
        if obj.ShowWaitingBar==1;
            close(hWaitingBar)
        end
    else
        % No separation is needed. Do prediction directly
        
        % Normalization of predictions points if needed
        if obj.NormInput==1
            for iInput = 1:obj.nInputVar
                input(:,iInput) = obj.scale(input(:,iInput), min(obj.InputData_True(:,iInput)),max(obj.InputData_True(:,iInput)),0,1);
            end
        end
        
        nInput = size(input,1);
        
        % Atual Prediction
        doPredictionUniversalKriging;
        
        % Save the krigingWeights
        obj.Weights = krigingWeights;
%         CovVecTot = obj.CovargramVectors;
        OutputTotal = Output;
    end
    
    % Save Covariogram vector (Needed for full facatorial and central
    % composite design) 
%     obj.CovargramVectors = CovVecTot;
    
%     if obj.NormOutput
%         % Normalization a = b*(max(out)-min(out))/(1-0) + min(out)
%         OutputTotal(:,1) = obj.scale(OutputTotal(:,1),0,1,...
%                        min(obj.getOutputData),max(obj.getOutputData));
%         % Do not scale standard deviation as for the calculation no
%         % normalized output values are used!!!!!!!!!!!!!!!!!!
%     end
    if obj.NormOutput
        % Scale = (x-minOld)*((maxNew-minNew)/(maxOld-minOld)) + minNew
        %       = (x-a)*b + c
        % -> x~N((y-a)*b+c,sigma^2*b^2)
        OutputTotal(:,1) = obj.scale(OutputTotal(:,1),0,1,...
                       min(obj.getOutputData),max(obj.getOutputData));
        b = (max(obj.getOutputData)-min(obj.getOutputData))/(1-0);
        OutputTotal(:,2) = OutputTotal(:,2)*b;
    end

%% --------------------- Nested Functions ---------------------------------
% -------------------------------------------------------------------------
    function []=doPredictionUniversalKriging()
        % Initialize Covariance Matrix
        obj.CovargramVectors = zeros(obj.nExperiments+obj.nBasisFct,nInput);
        
        % Calculate the distances between the new points and the provides
        % measurements and its corresponding covariance
        [~,index] = ndgrid(1:size(obj.InputData(:,:),1),1:nInput);
        index=index(:);
        if obj.CovariogramUsesEuclideanDistance
            obj.CovargramVectors=calcCovargramVectorEuclideanDistance(index);
        elseif obj.CovariogramUsesAbsoluteDistance
            obj.CovargramVectors=calcCovargramVectorAbsoluteDistance(index);
        else
            error('EitherCovariogramUsesEuclideanDistance or CovariogramUsesAbsoluteDistance have to be set true')
        end
        
        % Extend by basis function evaluiations
        for iBasis1 = 1 : obj.nBasisFct
            basis = obj.BasisFct{iBasis1}(obj.BasisFctParameters,nonNormInput(:,:));
            obj.CovargramVectors(obj.nExperiments+iBasis1,1:nInput) = (basis)';
        end

        % Calculate the extended variogram matrix and inverse of the
        % extended covariogram matrix, if not yet done
        if (obj.checkVariogram==0)
            obj.calcCovariogramMatrix();
            obj.checkVariogram = 1;
        end
        
        % ------------------------------------------------------------------------
        % Calculate Actual Prediction
        if obj.UseInverse
            % Inverse is direclty used
            
            % After explicite determination of the inverse of the extended
            % covariogram matrix
            if isempty(obj.InvCovariogramMatrix)
                obj.calcCovariogramMatrix;
            end
            
            % Special Case: Simple Kriging: Constant is given and is not
            % automatically estimated
            if obj.nBasisFctParameters == 0&&obj.getUseSimpleKriging==1
                krigingWeights = obj.InvCovariogramMatrix(1:obj.nExperiments,1:obj.nExperiments)*obj.CovargramVectors(1:obj.nExperiments,:);
            else
                krigingWeights = obj.InvCovariogramMatrix*obj.CovargramVectors;
            end
            
            
        else
            % Instrad of using Inverse directly, use Gaussian elemination
            
            if obj.nBasisFctParameters == 0&&obj.getUseSimpleKriging==1
                krigingWeights = obj.CovariogramMatrix(1:obj.nExperiments,1:obj.nExperiments)\obj.CovargramVectors(1:obj.nExperiments,:);
            else
                % Use gaussian elemination instead of the inverting the matrix
                krigingWeights = obj.CovariogramMatrix\obj.CovargramVectors;
            end
        end
        
        % Final Prediction
        prediction = [obj.OutputData;zeros(size(krigingWeights,1)-obj.nExperiments,1)]'*krigingWeights;
        prediction = prediction';
        
        % In case of simple Kriging, add basis function evaluations
        if obj.nBasisFctParameters == 0&&obj.getUseSimpleKriging==1
            if obj.nBasisFct~=1
                error('In case of simpel Kriging (nBasisFctParameters=0), exact one basis function has to tbe defined. But %i functions are defined',obj.nBasisFct)
            end
            % If prediction is far away from sample points 
            % -->(1-krigingWeights(1:obj.nExperiments,:)'*ones(obj.nExperiments,1))~=1
            % -->Influence of basis function is increasing
            prediction = prediction + (1-krigingWeights(1:obj.nExperiments,:)'*ones(obj.nExperiments,1)).*obj.BasisFct{1}(obj.BasisFctParameters,ones(nInput,1));
            % !!!To modify for general basis functions. 
            % !!!Idea: Create Basis function evaluation vector "B" (nInputXnBasisFct).
            % !!!Each entry (i,j) is the evaluation at inputData(i,:) of basis
            % !!!function BasisFct{j}. The final prediction is then:
            % !!!prediction = prediction+(1-krigingWeights(1:obj.nExperiments,:)'*ones(obj.nExperiments,1)).*sum(B,2));
        end
        
        % ------------------------------------------------------------------------
        % Calculate the estimated prediction error at each Point
        sigmaEstimation=calcSigmaEstimation();
        
        % Final Output
        Output = [prediction,sigmaEstimation];
    end
% -------------------------------------------------------------------------
    function [CovargramVector] = calcCovargramVectorEuclideanDistance(index)
        % Calculate euclidean distance between points. index contains
        % indices(nCombinations,2) of pairwise compared points
        
        distance = repmat(obj.InputData(:,:),nInput,1)-input(index,:);
            % Euclidean norm
        distance = sqrt(sum(distance.^2,2));
            % Convert to the correct shape
        distance = reshape(distance,obj.nExperiments,nInput);
            % Use calculated distances for the determination of the
            % covariogram vector
        CovargramVector(1:obj.nExperiments,:) = obj.CovarModel(distance,obj.AbsoluteInterpolation);
    end
% -------------------------------------------------------------------------
    function [CovargramVector] = calcCovargramVectorAbsoluteDistance(index)
        % Calculate euclidean distance between points. index contains
        % indices(nCombinations,2) of pairwise compared points
        
        
        distance = repmat(obj.InputData(:,:),nInput,1)-input(index,:);
            % Absolute Value
        distance = abs(distance);
            % Convert to the correct shape
        d = ones(obj.nExperiments,nInput);
        for iN = 1:nInput
            for iV = 1:obj.nInputVar
               d(:,(iN-1)*obj.nInputVar+iV)=distance(obj.nExperiments*(iN-1)+1:obj.nExperiments*(iN),iV);
            end
        end
            % "distance" size: 
            % (nExperiments)X(nInputVar)X(nInput)
        distance = reshape(d,obj.nExperiments,obj.nInputVar,nInput);
        
        % Calculate variogram vector
        CovargramVector(1:obj.nExperiments,:) = reshape(obj.CovarModel(distance,obj.AbsoluteInterpolation),obj.nExperiments,nInput);
    end

% -------------------------------------------------------------------------
    function [sigmaEstimation]=calcSigmaEstimation()
        % Higher correlation coefficient when you consider exact the same
        % sample point -> usually higher covariance expressed by
        % obj.sigmaError
%         covariance_zero  = ones(size(obj.CovargramVectors,2),size(krigingWeights,2))*obj.CovarModel(zeros(1,size(obj.DistInput,2)),1);
        covariance_zero  = ones(size(obj.CovargramVectors,2),size(krigingWeights,2))*obj.CovarModel(zeros(1,size(obj.DistInput,2)),0);
        
        % This formula is derive based on the derivatives Lagrange
        % Multiplier function with respect to the lagrange multiplier and
        % the weights (=sum(weight_i*(var(x_i-x_j)))+sum(mu_j*basisfct_j(x_i))=var(x^hat-x_i))
        % ... Implementation note: Use size(krigingWeights,1) since its
        % length ahs already changed in the case of simple Kriging

        % Following equations are the same as:
        % 1. 
%         sigmaEstimation = covariance_zero-krigingWeights'*obj.CovargramVectors(1:size(krigingWeights,1),:);
        % 2.
%         sigmaEstimation = covariance_zero - 
%                           - 2*krigingWeights'*obj.CovargramVectors(1:size(krigingWeights,1),:)
%                           + krigingWeights'*obj.CovariogramMatrix(1:size(krigingWeights,1),1:size(krigingWeights,1))*krigingWeights

        if isempty(obj.KriKitObjNoise)||isempty(obj.KriKitObjNoise.getOutputData) % Check if an initialized noise model exist
            covariance_zero = covariance_zero + obj.sigmaError^2;
        else
            predLogNoise = obj.KriKitObjNoise.prediction(nonNormInput(:,:));
            covariance_zero = bsxfun(@plus,covariance_zero, exp(predLogNoise(:,1)) );
        end
        
        sigmaEstimation = covariance_zero - 2*krigingWeights'*obj.CovargramVectors(1:size(krigingWeights,1),:) +...
                              krigingWeights'*obj.CovariogramMatrix(1:size(krigingWeights,1),1:size(krigingWeights,1))*krigingWeights;
        
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
