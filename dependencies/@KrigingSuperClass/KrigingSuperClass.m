classdef KrigingSuperClass<handle
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.
    %KRIGINGSUPERCLASS Summary of this class goes here
    %   Detailed explanation goes here
    
    %% Public Members
    properties(GetAccess='public',SetAccess='public')
    end
    
    %% Private Members
    properties(GetAccess='private',SetAccess='private')

    end
    
    %% Protected Members
    properties(GetAccess='protected',SetAccess='protected')
        %% General Variables
        % Contains the Provided Parameter Set. If Normailization is
        % activated this matrix conains double between 0 and 1
        InputData  = [];
        % Contains the actual Provided Parameter Set (not effected by
        % normalization). 
        InputData_True  = [];
        % Contains the provided Objective Values. If Normailization is
        % activated this matrix conains double between 0 and 1
        OutputData = [];
        % Contains the actual provided Objective Values(not effected by
        % normalization). 
        OutputData_True = [];
        % CovariogramModelChoice 1: parameters=[theta,sigma] covariance = obj.sigma^2*exp(-1/2*distance.^2/obj.theta^2)
        % CovariogramModelChoice 2: parameters=[theta,sigma,sigmaError] covariance =sigmaError^2 + obj.sigma^2*exp(-1/2*distance.^2/obj.theta^2);
        % CovariogramModelChoice 3: parameters=[theta,p,sigma,sigmaError] covariance = sigmaError^2 + obj.sigma^2*exp(-1/2*sum(distance^obj.p),obj.theta.^2)));
        % CovariogramModelChoice 4: parameters=[theta,sigma,sigmaError] covariance = sigmaError^2 + obj.sigma^2*exp(-1/2*sum(distance^2),obj.theta.^2)));
        % CovariogramModelChoice 5: parameters=[theta,sigma,sigmaError] covariance = sigmaError^2 + obj.sigma^2*(1 + sqrt(3)*weightedEuclidean).*exp((-sqrt(3)*weightedEuclidean)), 
        %                           weightedEuclidean=sqrt(distance.^2/obj.theta^2);
        % CovariogramModelChoice 6: as  CovariogramModelChoice 5 but with
        %                           weightedEuclidean = sqrt(sum(distance^2/obj.theta.^2))));
        % CovariogramModelChoice 7: covariance = sigmaError^2 + obj.sigma^2*(1 + sqrt(5)*weightedEuclidean+5/3*weightedEuclidean.^2).*exp((- sqrt(5)*weightedEuclidean));
        %                           weightedEuclidean=sqrt(distance.^2/obj.theta^2);
        % CovariogramModelChoice 8: as  CovariogramModelChoice 7 but with
        %                           weightedEuclidean = sqrt(sum(distance^2/obj.theta.^2))));
        % sigmaError is only added if distance is a vector of zeros.
        % CovariogramModelChoice has to be set manually after setting
        % input and output data. No default value is set.
        CovariogramModelChoice = -1;
        % This variable is adjusted automatically calling setCovariogramModelChoice
        CovariogramUsesEuclideanDistance = false;
        % This variable is adjusted automatically calling setCovariogramModelChoice
        CovariogramUsesAbsoluteDistance = false;
        % Calculate the Covariogram Matrix if it doesn't exist already
        checkVariogram   = 0;
        % This variable is used for saving the covariance function value
        % for two identical points
        covariance_zero = [];
        % Use SimpleKriging (Just using a constant as trend function but to
        % not estimate a linear coefficient infront of it). By default =
        % UseSimpleKriging = 0;
        UseSimpleKriging = false;
        %% Flags
        
        % Decide if the inverse of covariance should be calculated or if
        % the kriging interpolation should done using Gaussian elemination
        UseInverse = false;
        
        % Decide which solver for the optimization problems
        % 1 ... fmincon (SQP)
        % 2 ... Open Source: Genetic Algorithm
        %       Can be downloaded from 
        %       http://www.shef.ac.uk/acse/research/ecrg/gat.html and
        %       http://codem.group.shef.ac.uk/index.php/ga-toolbox 
        %       Introduction Paper:
        %       Chipperfield, A. J., and P. J. Fleming. "The MATLAB genetic
        %       algorithm toolbox." Applied Control Techniques Using
        %       MATLAB, IEE Colloquium on. IET, 1995.  
        % 3 ... GA (Genetic Algorithm) - default
        UseSolver = 3;
        % Option of the solver used for the optimization. Please consider
        % the help function for the specific solver to make sure that
        % "OptimizationOption" is of the correct form., 
        % E.g. in the case of fmincon and in Matlab2014 use 
        % - OptimizationOption = optimoptions(@fmincon,...)
        % But in Matlab2012 also in the case of fmincon use
        % - OptimizationOption = optimset(...)
        OptimizationOption = {};
        % CovariogramEstimationType = 1 ... Least Square Fitting of the
        %                                 variogram is performed (default)
        % CovariogramEstimationType = 2 ... CrossValidation is performed
        % CovariogramEstimationType = 3 ... Marginal Likelihood
        CovariogramEstimationType = 3;
        % Decide if the Input data shall be normalized (1,2) or not (0).
        % 0 ... Orginal InputData are used
        % 1 ... Normalization via to the range [0,1]
        % 2 ... Gauss transformation of the input data
        % By Default true
        NormInput = true;
        % Decide if the output data shall be normalized (1,2) or not (0).
        % 0 ... Orginal output data are used
        % 1 ... Normalization via to the range [0,1]
        % 2 ... Gauss transformation of the output data
        % By Default false
        NormOutput = false;
        % If AbsoluteInterpolation = 1 then the the prediction will be the
        % mean values of the samples at the measurment points. This may
        % lead to non-continous curves. By default = 0
        AbsoluteInterpolation = 0;
        % Prediction gets very slow when the input matrix is getting too
        % big. MaxSizeOfPredictions is the maximal number of points at
        % which the prediction is parallel calculated. By default 100
        MaxSizeOfPredictions = 100;
        % Decide if the waitbar should be showed. E.g. during prediction
        ShowWaitingBar = 1;
        %% Numbers of Variables
        % Number of Outputs that are provided for the estimation of the covariance matrix
        nExperiments    = 0;
        % Number of Inputs that are provided for the estimation of the covariance matrix
        nInputVar       = 0;
        % Number of Outputs that are provided for the estimation of the covariance matrix
        nOutputVar      = 0;
        %% CovarioGramm
        % Matrix containing the distances between the sample points (nInput x nInput)
        DistInput        = [];
        % Matrix containing the estimation of the variance values between
        % the inputs a function of the distance between them.
        VarianceEstimation = [];
        % Matrix containing the estimation of the correlation coefficient 
        % between the inputs a function of the distance between them
        CorrelationMatrix = [];
        % The explicit inverse of the extended variogram matrix. Only needed if UseInverse = 1
        InvVariogram    = [];
        % The explicit inverse of the extended covariogram matrix. Only needed if UseInverse = 1
        InvCovariogramMatrix
        % The extended variogram matrix. Extended by the trend functions.
        Variogram       = [];
        % The extended covariogram matrix. Extended by the trend functions.
        CovariogramMatrix
        % The explicit inverse of the extended covariogram matrix. Only needed if UseInverse = 1
        InvCoVariogram  = [];
        % The extended covariogram matrix. Extended by the trend functions.
        CoVariogram     = [];

        %--------------------------------------
        %% Radial function parameters
        LBCovariogramModelParameters = [0,0,0];
        UBCovariogramModelParameters = [];
        % Exponential radial function = sigmaError + sigma*exp(-(x/theta)^p)
        sigma=1;
        % Exponential radial function = sigmaError + sigma*exp(-(x/theta)^p)
        theta=1;
        % Exponential radial function = sigmaError + sigma*exp(-(x/theta)^p)
        p=2;
        % Exponential radial function = sigmaError + sigma*exp(-(x/theta)^p)
        sigmaError = 1e-10;
        % Squared difference between model of variance and the
        % actual empirical estimation
        ModelErrorVar = inf;
        % Initial value of the parameter of the variogram
        InitialCovariogramParameters = [];
        %% Universal Kriging
        % The number of basis functions saved in BasisFct
        nBasisFct = 0;
        % Degree of polynomial in bsis function
        MaxDegree = [];
        % A cell strcuture which saves the diffferent basis functions
        BasisFct = {};
        % The number of parameters in your basis functions, e.g. for
        % f=@(a,b,x)=a*x(1)+b has BasisFctParameters=2
        nBasisFctParameters = 0;
        % Number of Parameter for the current chosen covariogram
        nCovariogramParameters = 0;
        % Describe of type the basis function are, e.g. a*sin(b*x(1))+c has
        % BasisFctType='sinus'. 
        BasisFctType ='';
        % Different types of basis function which are already implemented
        % User means that a user defined function is applied
        allowed = {'polynomial','user'};
        
        % Array which saves the parameters of the defined basis function.
        % The particular use of the parameters depends on the function
        % itself. For more details regarding the method "LeastSquareFunctionBasisFct"
        BasisFctParameters = [];
        % Array which saves the coefficients c0-cn of the defined basis
        % functions.:(c0+c1*BasisFct_1+...+cn*BasisFct_n)
        BasisFctCoefficients = [];
        % Squared difference between basis function output and the
        % actual measurements
        ModelErrorBasisFct = inf;
        % Initial value of the parameters in the defined basis
        % function.
        InitialBasisFctParameters = [];
        % Lower bounds of the paramerters of the basis function
        LBBasisFctParameters = [];
        % Upper bounds of the paramerters of the basis function
        UBBasisFctParameters = [];
        % Decide if initial population should be uniform distributed over
        % the entire defined parameter space. LBBasisFctParameters and
        % UBBasisFctParameters must be defined here
        UseUniformDistributedPopulation = 0;
        %
        UseParallel = 1;
        %% Genetic Algorithm
%                 CrossoverFraction = 0.8
        
        % GA-Alg.: The fraction of the population at the next
        % generation, not including elite children, that is created by the
        % crossover function. By default 0.8
        CR = 0.8;
        % GA/DE-Alg.: Positive integer specifying the maximum number of
        % iterations before the algorithm holds
        % By default 100
        Generations = 100;
        % GA/DE-Alg.: Maximum of Time (in seconds) which pass until the optimization algorithm stops in seconds
        % By default 60
        TimeLimit = 60;
        % GA-Alg.: Size of the population. By default min(max(10*nInputVar,40),100)
        PopulationSize = -1;
        % GA-Alg.: Decide if the details shall be shown (=1) or
        % not (=0). % By default 0
        ShowDetails = 0;
        % GA-Alg.: Decide in which range the genes of the initial
        % population is allowed to be. PopInitMatrix is a matrix of size
        % nParameters x 2
        % The first column contains the maximum and the 2nd the minimum
        % values. If the matrix does not match the condition than the
        % setting from matlab are used. Does only work for the genetic
        % algorithm impletemeted in the matlab optimization toolbox and is
        % ignored if lower or upper bounds for the parameter are set
        PopInitMatrix = [];
        % Maximum number of Iterations of the local numerical solver such
        % as fmincon and fminsearch. By default nIterationsSolver = 1e3
        nIterationsSolver = 1e3;
        %% Cross Validation
        % Cross-Validation: Decide which of the quality calculation shall be used for the
        % calculation
        % 1 :(1-1/nExp*sum(((output_true-output_predicted)/sigma_predicted)^2))^2
        % 2 : sum((output_true-output_predicted)/^2)
        % 3 : sum(((output_true-output_predicted)/sigma_predicted)^2
        % 4 : Combines option 2 and 3 by  qualityTotal = quality1 + quality2/nUniqueInput/max(OutputBackup)^2;
        %     ... modification of quality2 makes sure that both qualtiy
        %     components are in the same numerical scale
        % By default 4
        CrossValidationType = 4;
        %% Ordinary regression Kriging
            % Kriging Object which is used when Kriging interpolation is calculated
        OrdinaryRegressionKrigingObject = [];
            % Interpolation using ordinary kriging at the measurement
            % data (with nugget)
        EstimatedOutputDataViaOrdinaryKriging = [];
        %% Additional information
        % Save Weights of Kriging Interpolation the last m entries contain
        % the lagrange multiplier for the m-basis functions
        Weights = [];
        % Save the covariagram vector used for the interpolation
        CovargramVectors = [];
        % Save the covariagram matrix of the predictions
        CovarMatrixOfPredictions = [];
        % Local Parameter Estimation
        %% Gaussian Process Regression 
        % Decide if Matlab GPR should be use (needs Matlab version >2015b and Statistic Toolbox)
        UseMatlabRegressionGP = false;
        GPR_Model = [];
        
    end
    
    methods
        %% Constructor
        function obj = KrigingSuperClass()
        end
        %% General Methods
        % ----------------------------------------------------------------
        function copy = copyObj(obj)
        % Create a shallow copy of the calling object.
            copy = eval(class(obj));
            meta = eval(['?',class(obj)]);
            for p = 1: size(meta.Properties,1)
                    pname = meta.Properties{p}.Name;
                try
                    eval(['copy.',pname,' = obj.',pname,';']);
                catch
                    error(['\nCould not copy ',pname,'.\n']);
%                     fprintf(['\nCould not copy ',pname,'.\n']);
                end
            end
        end
        %--------------------------------------
        % Define Basis functions
        setBasisFct(obj,type,varargin)
        %--------------------------------------
        % Calculate the covariance using the chosen model
        covariance = CovarModel(obj,distance,varargin);
        %--------------------------------------
        % Calculate the euclidean/absolute distance between the Inputs and
        % estimate the covariance function of the distance between inputs
        estimateVariance(obj)
        %--------------------------------------
        % Calculate (extended) Variogram Matrix
        calcCovariogramMatrix(obj);
        %--------------------------------------
        % Solve Least Square Problem
        SQ_CoVar = LeastSquareFunctionCovariogram(obj,varargin);
        % Solve Least Square Problem of finding the basis function
        % parameters
        SQ_CoVar = LeastSquareFunctionBasisFct(obj,varargin);
        % Maximization of marginal likelihood
        SQ_CoVar = MarginalLikelihoodCovarPara(obj,varargin);
        % Solve Least Square Problem: fitting the basis function parameter (Local via fmincon/fminSearch)
        solveLeastSquareBasisFct(obj);
        % Solve Least Square Problem: fitting the basis function parameter (Global via Genetic Algorithm)
        solveLeastSquareBasisFctGA(obj);
        % Solve Least Square Problem: fitting the basis function parameter (Global via OpenCscource Algorithm)
        solveLeastSquareBasisFctGA2(obj);
        % Solve Least Square Problem of the CoVar model (Local via fmincon)
        solveLeastSquareCovariogram(obj);
        % Solve Least Square Problem of the CoVar model (Global via Genetic Algorithm - Optimization Toolbox)
        solveLeastSquareCovariogramGA(obj);
        % Solve Least Square Problem of the CoVar model (Global via Genetic Algorithm Open Scource)
        solveLeastSquareCovariogramGA2(obj);
        % Generate RGP model for further Kriging estimation
        generateRegressionGPModel(obj);
        % gaCreationUniformFeasible
        function[Population]=gaCreationUniformFeasible(obj,GenomeLength,FitnessFcn,options)
            Population = bsxfun(@times,rand(sum(options.PopulationSize),GenomeLength),(obj.getUBBasisFctParameters - obj.getLBBasisFctParameters));
            Population = bsxfun(@plus,Population,obj.getLBBasisFctParameters);
        end
        %--------------------------------------
        % This function makes all the work which have to be done before you
        % can perform the prediction
        %--------------------------------------
        []=makePrework(obj);
        % This function calculates the quality covariogram parameters. The
        % quality is calculated as
        [quality]=calcCrossOverQuality(obj,varargin);
        % -----------------------------------------------------------------
        % Prediction of the output for the given input points using the
        % Covariogram matrix
        [output]=prediction(obj,input);
        %--------------------------------------
        % Plot Variogramm
        plotVariogram(obj);
        % ----------------------------------------------------------------
        function f=scale(obj,x,minX,maxX,LB,UB)
            % f=scale(x,minX,maxX,LB,UB)
            % Scales x which is a variable in the range [minX,maxX] to the 
            % range [LB,UB]
            
            if (maxX-minX)<1e-10
                error('During Sacling: maxX must be bigger than minX. Not Smaller and also not equal!')
            end
            f=(x-minX)./(maxX-minX).*(UB-LB)+LB;
        end
        % ----------------------------------------------------------------
        estimateBasisFctCoefficients(obj,varargin)
        % ----------------------------------------------------------------
        [] = estimateBasisParametersViaOrdinaryRegressionKriging(obj)
        % ----------------------------------------------------------------
        [] = plotResultsOfOrdinaryRegressionKriging(obj)
        % ----------------------------------------------------------------
        [options]=setOptionsGA(obj)
        % ----------------------------------------------------------------
        [options]=setOptionsFminCon(obj)
        % ----------------------------------------------------------------
        [options]=setOptionsFminUnc(obj)
        % ----------------------------------------------------------------
        [options]=setOptionsFminSearch(obj)
        % ----------------------------------------------------------------
        [] = InitializeCovariogramParameters(obj)
        % ----------------------------------------------------------------
        [] = saveOptimizedVariogramParameters(obj,optPara)
        % ----------------------------------------------------------------
%         [] = defineParametersInLSQFunction (obj,varargin)
        % ----------------------------------------------------------------
        [realizationCurve,nonNormSampleLocations] = doConditionalSimulation_SeqSim(obj,varargin)
        % ----------------------------------------------------------------
        [realizationCurve,nonNormSampleLocations] = doConditionalSimulation_ResAlgInEq(obj,varargin)
        % ----------------------------------------------------------------
        function [] = checkIfNumberOfParameterValuesIsCorrect(obj,varargin)
            if(length(varargin{1})>obj.nCovariogramParameters)
                error('Covariogram needs only %i parameter but %i are defined!',obj.nCovariogramParameters,length(varargin{1}));
            end
            if(length(varargin{1})<obj.nCovariogramParameters)
                error('Covariogram needs %i parameter but only %i are defined!',obj.nCovariogramParameters,length(varargin{1}));
            end
        end
        %% Get Functions
        function [DistInput]=getDistInput(obj) 
            DistInput = obj.DistInput;
        end
        %--------------------------------------
        function [VarianceEstimation]=getVarianceEstimation(obj) 
            VarianceEstimation = obj.VarianceEstimation;
        end
        %--------------------------------------
        function [nInputVar]=getnInputVar(obj) 
            nInputVar = obj.nInputVar;
        end
        %--------------------------------------
        function [nExperiments]=getnExperiments(obj) 
            nExperiments = obj.nExperiments;
        end
        %--------------------------------------
        function [parameters]=getCovariogramModelParameters(obj) 
            % [parameters]=getCovariogramModelParameters()
            % This function gives out the parameters of the model which is
            % used for the covariance
            % CovariogramModelChoice 1: parameters=[theta,sigma] covariance = obj.sigma^2*exp(-1/2*distance.^2/obj.theta^2)
            % CovariogramModelChoice 2: parameters=[theta,sigma,sigmaError] covariance =sigmaError + obj.sigma^2*exp(-1/2*distance.^2/obj.theta^2);
            % CovariogramModelChoice 3: parameters=[theta,p,sigma,sigmaError] covariance = sigmaError + obj.sigma^2*exp(-1/2*sum(distance^obj.p),obj.theta.^2)));
            % CovariogramModelChoice 4: parameters=[theta,sigma,sigmaError] covariance = sigmaError + obj.sigma^2*exp(-1/2*sum(distance^2),obj.theta.^2)));
            % CovariogramModelChoice 5: parameters=[theta,sigma,sigmaError] covariance = sigmaError + obj.sigma^2*(1 + sqrt(3)*weightedEuclidean).*exp((-sqrt(3)*weightedEuclidean)), 
            %                           weightedEuclidean=sqrt(distance.^2/obj.theta^2);
            % CovariogramModelChoice 6: as  CovariogramModelChoice 5 but with
            %                           weightecalcInterpolation_3DdEuclidean = sqrt(sum(distance^2/obj.theta.^2))));
            % CovariogramModelChoice 7: covariance = sigmaError + obj.sigma^2*(1 + sqrt(5)*weightedEuclidean+5/3*weightedEuclidean.^2).*exp((- sqrt(5)*weightedEuclidean));
            %                           weightedEuclidean=sqrt(distance.^2/obj.theta^2);
            % CovariogramModelChoice 8: as  CovariogramModelChoice 7 but with
            %                           weightedEuclidean = sqrt(sum(distance^2/obj.theta.^2))));
            switch obj.CovariogramModelChoice
                case 1
                    parameters = [obj.theta,obj.sigma];
                case {2,5,7}
%                     parameters = [obj.theta,obj.sigma,obj.p,obj.sigmaError];
                    parameters = [obj.theta,obj.sigma,obj.sigmaError];
                case {3}
                    parameters = [obj.theta,obj.p,obj.sigma,obj.sigmaError];
                case {4,6,8}
                    parameters = [obj.theta,obj.sigma,obj.sigmaError];
                otherwise
                    error('CovariogramModelChoice=%i is not defined',obj.CovariogramModelChoice);
            end
        end
        %--------------------------------------
        function [ModelErrorVar] = getModelErrorVar(obj)
            ModelErrorVar=obj.ModelErrorVar;
        end
        %--------------------------------------
        function [ModelErrorBasisFct] = getModelErrorBasisFct(obj)
            ModelErrorBasisFct=obj.ModelErrorBasisFct;
        end
        %--------------------------------------
        function [Variogram] = getVariogramMatrix(obj)
            Variogram = obj.Variogram;
        end
        %--------------------------------------
        function [CovariogramMatrix] = getCovariogramMatrix(obj)
            CovariogramMatrix = obj.CovariogramMatrix;
        end
        %--------------------------------------
        function [CovariogramModelChoice] = getCovariogramModelChoice(obj)
            CovariogramModelChoice=obj.CovariogramModelChoice;
        end
        %--------------------------------------
        function [Input] = getInputData(obj)
            Input = obj.InputData_True;
        end
        %------------------------------------------------------------------
        function [Output] = getOutputData(obj)
            Output = obj.OutputData_True;
        end
        %------------------------------------------------------------------
        function [LowerBound] = getLBCovariogramModelParameters(obj)
            % CovariogramModelChoice 1: parameters=[theta,sigma] covariance = obj.sigma^2*exp(-1/2*distance.^2/obj.theta^2)
            % CovariogramModelChoice 2: parameters=[theta,sigma,sigmaError] covariance =sigmaError + obj.sigma^2*exp(-1/2*distance.^2/obj.theta^2);
            % CovariogramModelChoice 3: parameters=[theta,p,sigma,sigmaError] covariance = sigmaError + obj.sigma^2*exp(-1/2*sum(distance^obj.p),obj.theta.^2)));
            % CovariogramModelChoice 4: parameters=[theta,sigma,sigmaError] covariance = sigmaError + obj.sigma^2*exp(-1/2*sum(distance^2),obj.theta.^2)));
            % CovariogramModelChoice 5: parameters=[theta,sigma,sigmaError] covariance = sigmaError + obj.sigma^2*(1 + sqrt(3)*weightedEuclidean).*exp((-sqrt(3)*weightedEuclidean)), 
            %                           weightedEuclidean=sqrt(distance.^2/obj.theta^2);
            % CovariogramModelChoice 6: as  CovariogramModelChoice 5 but with
            %                           weightedEuclidean = sqrt(sum(distance^2/obj.theta.^2))));
            % CovariogramModelChoice 7: covariance = sigmaError + obj.sigma^2*(1 + sqrt(5)*weightedEuclidean+5/3*weightedEuclidean.^2).*exp((- sqrt(5)*weightedEuclidean));
            %                           weightedEuclidean=sqrt(distance.^2/obj.theta^2);
            % CovariogramModelChoice 8: as  CovariogramModelChoice 7 but with
            %                           weightedEuclidean = sqrt(sum(distance^2/obj.theta.^2))));
            LowerBound = obj.LBCovariogramModelParameters;
        end
        %------------------------------------------------------------------
        function [UpperBound] = getUBCovariogramModelParameters(obj)
            UpperBound = obj.UBCovariogramModelParameters;
        end
        %------------------------------------------------------------------
%         function [CoVariance]=getCoVar(obj)
%             CoVariance = obj.CoVar;
%         end
%         %------------------------------------------------------------------
%         function [InvCoVariance]=getInvCoVar(obj)
%             InvCoVariance = obj.InvCoVar;
%         end
        function [CoVariogram]=getCoVariogram(obj)
            CoVariogram = obj.CoVariogram;
        end
        %------------------------------------------------------------------
        function [InvCoVariogram]=getInvCoVariogram(obj)
            InvCoVariogram = obj.InvCoVariogram;
        end
        %------------------------------------------------------------------
        function [nBasisFct] = getnBasisFct(obj)
        % Gives out the number of basis functions
            nBasisFct = obj.nBasisFct;
        end
        %------------------------------------------------------------------
        function [MaxDegree] = getMaxDegree(obj)
        % Gives out the number of basis functions
            MaxDegree = obj.MaxDegree;
        end
        %------------------------------------------------------------------
        function [BasisFct] = getBasisFct(obj)
            BasisFct = obj.BasisFct;
        end
        %------------------------------------------------------------------
        function [CR]=getCR(obj)
            % Get the fraction of the population at the next generation, not including elite children, that is created by the crossover function
            CR = obj.CR;
        end
        %------------------------------------------------------------------
        function [Generations]=getGenerations(obj)
            % Get the Positive integer specifying the maximum number of iterations before the algorithm halts
            Generations = obj.Generations;
        end
        %------------------------------------------------------------------
        function [nIterationsSolver]=getnIterationsSolver(obj)
            nIterationsSolver = obj.nIterationsSolver;
        end
%         %------------------------------------------------------------------
%         function [MigrationFraction]=getMigrationFraction(obj)
%             % Get the Scalar between 0 and 1 specifying the fraction of individuals in each subpopulation that migrates to a different subpopulation
%             MigrationFraction = obj.MigrationFraction;
%         end
        %------------------------------------------------------------------
        function [NormInput]=getNormInput(obj)
            % (1) if the input is normalized. (0) if not
            NormInput = obj.NormInput;
        end
        %------------------------------------------------------------------
        function [NormOutput]=getNormOutput(obj)
            % (1) if the Output is normalized. (0) if not
            NormOutput = obj.NormOutput ;
        end
        %------------------------------------------------------------------
        function [checkVariogram]=getcheckVariogram(obj)
            % (1) if the Output is normalized. (0) if not
            checkVariogram = obj.checkVariogram ;
        end
        %------------------------------------------------------------------
        function [TimeLimit]=getTimeLimit(obj)
            % Gives out the maximal time until the optimization algorithm stops
            TimeLimit = obj.TimeLimit ;
        end
        %
        %------------------------------------------------------------------
        function [MeanCoeffBasisFct]=getMeanCoeffBasisFct(obj)
            % Gives out the mean values of the coefficients of the basis functions used in universal kriging and bayesian kriging
            MeanCoeffBasisFct = obj.MeanCoeffBasisFct ;
        end
        %------------------------------------------------------------------
        function [CoVarCoeffBasisFct]=getCoVarCoeffBasisFct(obj)
            % Gives out the covariance matrix of the coefficients of the basis functions used in universal kriging and bayesian kriging
            CoVarCoeffBasisFct = obj.CoVarCoeffBasisFct ;
        end
        %------------------------------------------------------------------
        function [AverageProductOfOutputs]=getAverageProductOfOutputs(obj)
            AverageProductOfOutputs = obj.AverageProductOfOutputs ;
        end
        %------------------------------------------------------------------
        function [ProdWeightedFunctions]=getProdWeightedFunctions(obj)
            ProdWeightedFunctions = obj.ProdWeightedFunctions;
        end
        %------------------------------------------------------------------
        function [SumSumOfWeightedBasisFcts]=getSumSumOfWeightedBasisFcts(obj)
            SumSumOfWeightedBasisFcts = obj.SumSumOfWeightedBasisFcts;
        end
        %------------------------------------------------------------------
        function [UseSolver]=getUseSolver(obj)
            UseSolver = obj.UseSolver;
        end
        %------------------------------------------------------------------
        function [CorrelationMatrix]=getCorrelationMatrix(obj)
            CorrelationMatrix = obj.CorrelationMatrix;
        end
        %------------------------------------------------------------------
        function [nBasisFctParameters]=getnBasisFctParameters(obj)
            nBasisFctParameters = obj.nBasisFctParameters;
        end
        %------------------------------------------------------------------
        function [BasisFctType]=getBasisFctType(obj)
            BasisFctType = obj.BasisFctType;
        end
        %------------------------------------------------------------------
        function [InitialBasisFctParameters]=getInitialBasisFctParameters(obj)
            InitialBasisFctParameters = obj.InitialBasisFctParameters;
        end
        %------------------------------------------------------------------
        function [BasisFctParameters]=getBasisFctParameters(obj)
            BasisFctParameters = obj.BasisFctParameters;
        end
        %------------------------------------------------------------------
        function [LBBasisFctParameters]=getLBBasisFctParameters(obj)
            LBBasisFctParameters = obj.LBBasisFctParameters;
        end
        %------------------------------------------------------------------
        function [UBBasisFctParameters]=getUBBasisFctParameters(obj)
            UBBasisFctParameters = obj.UBBasisFctParameters;
        end
        %
        %------------------------------------------------------------------
        function [PopulationSize]=getPopulationSize(obj)
            PopulationSize = obj.PopulationSize;
        end
        %------------------------------------------------------------------
        function [ShowDetails]=getShowDetails(obj)
            ShowDetails = obj.ShowDetails;
        end
        %------------------------------------------------------------------
        function [nIterationCrossValidation]=getnIterationCrossValidation(obj)
            nIterationCrossValidation = obj.nIterationCrossValidation;
        end
        %------------------------------------------------------------------
        function [nLowerUpperBoundSections]=getnLowerUpperBoundSections(obj)
            nLowerUpperBoundSections = obj.nLowerUpperBoundSections;
        end
        %------------------------------------------------------------------
        function [nTimeSections]=getnTimeSections(obj)
            nTimeSections = obj.nTimeSections;
        end
        %------------------------------------------------------------------
        function [nCrossOverSections]=getnCrossOverSections(obj)
            nCrossOverSections = obj.nCrossOverSections;
        end
        %------------------------------------------------------------------
        function [nCrossoverFractionSections]=getnCrossoverFractionSections(obj)
            nCrossoverFractionSections = obj.nCrossoverFractionSections;
        end
        %------------------------------------------------------------------
        function [nGenerationSections]=getnGenerationSections(obj)
            nGenerationSections = obj.nGenerationSections;
        end
        %------------------------------------------------------------------
        function [CrossValidationType]=getCrossValidationType(obj)
            CrossValidationType = obj.CrossValidationType;
        end
        %------------------------------------------------------------------
        function [optionMatrix]=getoptionMatrix(obj)
            optionMatrix = obj.optionMatrix;
        end
        %------------------------------------------------------------------
        function [stringVector]=getstringVector(obj)
            stringVector = obj.stringVector;
        end
        %------------------------------------------------------------------
        function [CovariogramEstimationType]=getCovariogramEstimationType(obj)
            CovariogramEstimationType = obj.CovariogramEstimationType;
        end 
        %------------------------------------------------------------------
        function [BasisFctCoefficients]=getBasisFctCoefficients(obj)
            BasisFctCoefficients = obj.BasisFctCoefficients;
        end
        %------------------------------------------------------------------
        function [InitialCovariogramParameters]=getInitialCovariogramParameters(obj)
            InitialCovariogramParameters = obj.InitialCovariogramParameters;
        end
        %------------------------------------------------------------------
        function [OrdinaryRegressionKrigingObject]=getOrdinaryRegressionKrigingObject(obj)
            OrdinaryRegressionKrigingObject = obj.OrdinaryRegressionKrigingObject;
        end
        %------------------------------------------------------------------
        function [EstimatedOutputDataViaOrdinaryKriging] = getEstimatedOutputDataViaOrdinaryKriging(obj)
            EstimatedOutputDataViaOrdinaryKriging = obj.EstimatedOutputDataViaOrdinaryKriging;
        end
        %------------------------------------------------------------------
        function [Weights] = getWeights(obj)
            Weights = obj.Weights;
        end
        %------------------------------------------------------------------
        function [AbsoluteInterpolation]=getAbsoluteInterpolation(obj)
            AbsoluteInterpolation = obj.AbsoluteInterpolation;
        end
        %------------------------------------------------------------------
        function [MaxSizeOfPredictions]=getMaxSizeOfPredictions(obj)
             MaxSizeOfPredictions = obj.MaxSizeOfPredictions;
        end
        %------------------------------------------------------------------
        function [RelativePartOfCrossValidation]=getRelativePartOfCrossValidation(obj)
            RelativePartOfCrossValidation = obj.RelativePartOfCrossValidation;
        end
        %------------------------------------------------------------------
        function [RelativePartOfConfidenceInterval]=getRelativePartOfConfidenceInterval(obj)
            RelativePartOfConfidenceInterval = obj.RelativePartOfConfidenceInterval;
        end
        %------------------------------------------------------------------
        function [ShowWaitingBar]=getShowWaitingBar(obj)
            ShowWaitingBar = obj.ShowWaitingBar;
        end
        function [InvCovariogramMatrix]=getInvCovariogramMatrix(obj)
            InvCovariogramMatrix = obj.InvCovariogramMatrix;
        end
        %------------------------------------------------------------------
        function [UseUniformDistributedPopulation]=getUseUniformDistributedPopulation(obj)
            UseUniformDistributedPopulation = obj.UseUniformDistributedPopulation;
        end
        %------------------------------------------------------------------
        function [CovargramVectors]=getCovargramVectors(obj)
            CovargramVectors = obj.CovargramVectors;
        end
        function [CovarMatrixOfPredictions]=getCovarMatrixOfPredictions(obj)
            CovarMatrixOfPredictions = obj.CovarMatrixOfPredictions;
        end
        %------------------------------------------------------------------
        function [CoeffLSQVargramCrossValCumLSQ]=getCoeffLSQVargramCrossValCumLSQ(obj)
            CoeffLSQVargramCrossValCumLSQ = obj.CoeffLSQVargramCrossValCumLSQ;
        end
        %------------------------------------------------------------------
        function [PopInitMatrix]=getPopInitMatrix(obj)
            PopInitMatrix = obj.PopInitMatrix;
        end
        %------------------------------------------------------------------
        function [UseSimpleKriging]=getUseSimpleKriging(obj)
            UseSimpleKriging = obj.UseSimpleKriging;
        end
        %------------------------------------------------------------------
        function [UseParallel]=getUseParallel(obj)
            UseParallel = obj.UseParallel;
        end
        %------------------------------------------------------------------
        function [UseInverse]=getUseInverse(obj)
            % (1) if the Output is normalized. (0) if not
            UseInverse = obj.UseInverse;
        end
        %------------------------------------------------------------------
        function [nCovariogramParameters]=getnCovariogramParameters(obj)
            nCovariogramParameters = obj.nCovariogramParameters;
        end
        %------------------------------------------------------------------
        function [UseMatlabRegressionGP]=getUseMatlabRegressionGP(obj)
            UseMatlabRegressionGP = obj.UseMatlabRegressionGP;
        end
        %------------------------------------------------------------------
        function [GPR_Model]=getGPR_Model(obj)
            GPR_Model = obj.GPR_Model;
        end
        %------------------------------------------------------------------
        function [CovariogramUsesEuclideanDistance]=getCovariogramUsesEuclideanDistance(obj)
            CovariogramUsesEuclideanDistance = obj.CovariogramUsesEuclideanDistance;
        end
        %------------------------------------------------------------------
        function [CovariogramUsesAbsoluteDistance]=getCovariogramUsesAbsoluteDistance(obj)
            CovariogramUsesAbsoluteDistance = obj.CovariogramUsesAbsoluteDistance;
        end

        %% Set Functions
        function [] = setInputData(obj,Input)
            % Define the InputVariables the number of rows should be equal
            % to the number of experiments and the columns equal to the
            % number of input variables. By setting the NormOutput flag the
            % InputData are getting normalized
            % NormInput=0 ... Orginal InputData are used
            % NormInput=1 ... Normalization via to the range [0,1]
            % NormInput=2 ... Gauss transformation of the input data
            % 
            % You can set: 
            % NormInput ... if true, input data are scaled between 0 and 1.
            %               This makes prediction more robust
            %
            
            if size(Input,1)==1&&size(Input,2)>=1
                Input = Input';
                warning('InputData seem to be a row vector and is transposed')
            end
            if any(any(isinf(Input)))||any(any(isnan(Input)))
                error('No "inf" or "nan" entry allowed in input data')
            end
            
            obj.InputData = Input;
            obj.InputData_True = Input;
            obj.nExperiments = size(Input,1);
            obj.nInputVar = size(Input,2);
            if(obj.nInputVar>obj.nExperiments)
                warning('InputData: More Input variables as Experiments');
            end
            
            
            
            
            
            switch obj.NormInput
                case 0
                case 1
                    % Normalization to the range [0,1]
                    for iInput = 1:obj.nInputVar
                        if min(obj.InputData(:,iInput))==max(obj.InputData(:,iInput))
                            error('InputVariable %i does never change its value',iInput)
                        end
                        obj.InputData(:,iInput) = obj.scale(obj.InputData(:,iInput),min(obj.InputData(:,iInput)),...
                            max(obj.InputData(:,iInput)),0,1);
                    end
                case 2
                    % Get emperic cumulativeFrequency which is distributed
                    % [0,1] but following an unknown distribution
                    [cumulativeFrequency]=ecdf(obj.InputData);
                    
                    % -0.5/n since otherwise the last probability would be
                    % 1 and therefor the inverser would be infinity
                    cumulativeFrequency = cumulativeFrequency - 0.5/length(obj.InputData);
                    
                    % The first component has the probability of zero
                    cumulativeFrequency = cumulativeFrequency(2:end);
                    obj.InputData = icdf('Normal',cumulativeFrequency,0,1);
                otherwise
                    error('NormInput = %i is not defined',obj.NormInput);
            end
        end
        %------------------------------------------------------------------
        % Define the outputVariables
        function [] = setOutputData(obj,Output)
            % Define the OutputVariables the number of rows should be equal
            % to the number of experiments and the columns equal to one By  
            % setting the NormOutput flag the output data are normalized
            % NormOutput=0 ... Orginal InputData are used
            % NormOutput=1 ... Normalization via to the range [0,1]
            % NormOutput=2 ... Gauss transformation of the Output data
            if(size(Output,2)~=1&&size(Output,1)~=1)
               error('Dimension of Output vector is not correct! The number of columns should be 1 but it is %i',size(Output,2));
            end
            if(size(Output,2)~=1&&size(Output,1)==1)
               Output=Output';
               warning('Output data is provided in a row vector, a column vector is expected. The provided vector is transposed.');
            end
            if(size(Output,1)~=obj.nExperiments)
               error('The number of entries in the InputData (%i) and in the OutputData (%i) are not the same! \nIn both vectors the rows should be associated with the particular experiment. Please define the input data first',obj.nExperiments,size(Output,1));
            end
            if any(any(isinf(Output)))||any(any(isnan(Output)))
                error('No "inf" or "nan" entry allowed in output data')
            end
            
            obj.OutputData = Output;
            obj.OutputData_True = Output;
            obj.nOutputVar = obj.nExperiments;
            
            switch obj.NormOutput
                case 0
                    
                case 1
                    obj.OutputData = obj.scale(obj.OutputData,...
                    min(obj.OutputData),max(obj.OutputData),0,1);
                case 2
                    % Get emperic cumulativeFrequency which is distributed
                    % [0,1] but following an unknown distribution
                    [cumulativeFrequency]=ecdf(obj.OutputData);
                    
                    % -0.5/n since otherwise the last probability would be
                    % 1 and therefor the inverser would be infinity
                    cumulativeFrequency = cumulativeFrequency - 0.5/length(obj.OutputData);
                    
                    % The first component has the probability of zero
                    cumulativeFrequency = cumulativeFrequency(2:end);
                    obj.OutputData =  icdf('Normal',cumulativeFrequency,0,1);
                otherwise
                    error('NormOutput = %i is not defined',obj.NormOutput);
            end
                    
%             if obj.NormOutput==1
%                     obj.OutputData = obj.scale(obj.OutputData,...
%                         min(obj.OutputData),max(obj.OutputData),0,1);
%             end
        end
        %------------------------------------------------------------------
        % Choose the vlaues for the exponetial radial function
        function []=setCovariogramModelParameters(obj,varargin) 
            % []=setCovariogramModelParameters(varargin) 
            % This function set the Kriging parameters. The Parameters
            % depends on your model choice for the covariance
            % CovariogramModelChoice 1: parameters=[theta,sigma] covariance = obj.sigma^2*exp(-1/2*distance.^2/obj.theta^2)
            % CovariogramModelChoice 2: parameters=[theta,sigma,sigmaError] covariance =sigmaError + obj.sigma^2*exp(-1/2*distance.^2/obj.theta^2);
            % CovariogramModelChoice 3: parameters=[theta,p,sigma,sigmaError] covariance = sigmaError + obj.sigma^2*exp(-1/2*sum(distance^obj.p),obj.theta.^2)));
            % CovariogramModelChoice 4: parameters=[theta,sigma,sigmaError] covariance = sigmaError + obj.sigma^2*exp(-1/2*sum(distance^2),obj.theta.^2)));
            % CovariogramModelChoice 5: parameters=[theta,sigma,sigmaError] covariance = sigmaError + obj.sigma^2*(1 + sqrt(3)*weightedEuclidean).*exp((-sqrt(3)*weightedEuclidean)), 
            %                           weightedEuclidean=sqrt(distance.^2/obj.theta^2);
            % CovariogramModelChoice 6: as  CovariogramModelChoice 5 but with
            %                           weightedEuclidean = sqrt(sum(distance^2/obj.theta.^2))));
            % CovariogramModelChoice 7: covariance = sigmaError + obj.sigma^2*(1 + sqrt(5)*weightedEuclidean+5/3*weightedEuclidean.^2).*exp((- sqrt(5)*weightedEuclidean));
            %                           weightedEuclidean=sqrt(distance.^2/obj.theta^2);
            % CovariogramModelChoice 8: as  CovariogramModelChoice 7 but with
            %                           weightedEuclidean = sqrt(sum(distance^2/obj.theta.^2))));

            
            obj.checkIfNumberOfParameterValuesIsCorrect(varargin{1})
            
            switch obj.CovariogramModelChoice
                case 1
%                     if(length(varargin{1})>2)
%                         warning('The exponential model needs only 2 parameter but %i are defined!',length(varargin{1}));
%                     end
%                     if(length(varargin{1})<2)
%                         error('The exponential model needs 2 parameter but only %i are defined!',length(varargin{1}));
%                     end
                    obj.theta       = varargin{1}(1);
                    obj.sigma       = varargin{1}(2);
                    
%                     obj.CovariogramUsesEuclideanDistance = true;
%                     obj.p           = varargin{1}(3);
                case {2,5,7}
                    obj.theta       = varargin{1}(1);
                    obj.sigma       = varargin{1}(2);
                    obj.sigmaError  = varargin{1}(3);
                case 3
                    nTotal = obj.nInputVar;
                    obj.theta       = varargin{1}(1:nTotal);
                    obj.p           = varargin{1}(nTotal+1:2*nTotal);
                    obj.sigma       = varargin{1}(2*nTotal+1);
                    obj.sigmaError  = varargin{1}(2*nTotal+2);
                case {4,6,8}
                    nTotal = obj.nInputVar;
                    obj.theta       = varargin{1}(1:nTotal);
                    obj.sigma       = varargin{1}(nTotal+1);
                    obj.sigmaError  = varargin{1}(nTotal+2);
                otherwise
                    error('Non-acceptable model choice. The parameter CovariogramModelChoice = %i is not allowed',obj.CovariogramModelChoice);
            end
            
%             switch obj.CovariogramModelChoice
%                 case {1,2,5,7}
%                     obj.CovariogramUsesEuclideanDistance = true;
%                     obj.CovariogramUsesAbsoluteDistance = false;
%                 case {3,4,6,8}
%                     obj.CovariogramUsesAbsoluteDistance = true;
%                     obj.CovariogramUsesEuclideanDistance = false;
%                 otherwise
%                     error('CovariogramModelChoice = %d is not defined',obj.CovariogramModelChoice)
%             end
            
            try
                obj.calcCovariogramMatrix();
            catch exception
                warning('Covariogram Matrix could be recalculated! Do it manually')
                % Rethrow original error.
                rethrow(exception)
            end
        end
        %------------------------------------------------------------------
        function [] = setLBCovariogramModelParameters(obj,varargin)
            % [] = setLBCovariogramModelParameters(varargin)
            % Notice this function shall only be called AFTER
            % setInput/setCovariogramModelChoice
            % CovariogramModelChoice 1: parameters=[theta,sigma] covariance = obj.sigma^2*exp(-1/2*distance.^2/obj.theta^2)
            % CovariogramModelChoice 2: parameters=[theta,sigma,sigmaError] covariance =sigmaError + obj.sigma^2*exp(-1/2*distance.^2/obj.theta^2);
            % CovariogramModelChoice 3: parameters=[theta,p,sigma,sigmaError] covariance = sigmaError + obj.sigma^2*exp(-1/2*sum(distance^obj.p),obj.theta.^2)));
            % CovariogramModelChoice 4: parameters=[theta,sigma,sigmaError] covariance = sigmaError + obj.sigma^2*exp(-1/2*sum(distance^2),obj.theta.^2)));
            % CovariogramModelChoice 5: parameters=[theta,sigma,sigmaError] covariance = sigmaError + obj.sigma^2*(1 + sqrt(3)*weightedEuclidean).*exp((-sqrt(3)*weightedEuclidean)), 
            %                           weightedEuclidean=sqrt(distance.^2/obj.theta^2);
            % CovariogramModelChoice 6: as  CovariogramModelChoice 5 but with
            %                           weightedEuclidean = sqrt(sum(distance^2/obj.theta.^2))));
            % CovariogramModelChoice 7: covariance = sigmaError + obj.sigma^2*(1 + sqrt(5)*weightedEuclidean+5/3*weightedEuclidean.^2).*exp((- sqrt(5)*weightedEuclidean));
            %                           weightedEuclidean=sqrt(distance.^2/obj.theta^2);
            % CovariogramModelChoice 8: as  CovariogramModelChoice 7 but with
            %                           weightedEuclidean = sqrt(sum(distance^2/obj.theta.^2))));
            obj.checkIfNumberOfParameterValuesIsCorrect(varargin{1})
            if size(varargin{1},1)>1&&size(varargin{1},2)==1
                varargin{1} = varargin{1}';
            end
            obj.LBCovariogramModelParameters = varargin{1};
        end
        %------------------------------------------------------------------
        function [] = setUBCovariogramModelParameters(obj,varargin)
            % [] = setUBCovariogramModelParameters(varargin)
            % Notice this function shall only be called AFTER
            % setInput/setCovariogramModelChoice
            
            obj.checkIfNumberOfParameterValuesIsCorrect(varargin{1})
            if size(varargin{1},1)>1&&size(varargin{1},2)==1
                varargin{1} = varargin{1}';
            end
            obj.UBCovariogramModelParameters = varargin{1};
        end
        %------------------------------------------------------------------
        function []=setCovariogramModelChoice(obj,varargin) 
            % CovariogramModelChoice 1: parameters=[theta,sigma] covariance = obj.sigma^2*exp(-1/2*distance.^2/obj.theta^2)
            % CovariogramModelChoice 2: parameters=[theta,sigma,sigmaError] covariance =sigmaError + obj.sigma^2*exp(-1/2*distance.^2/obj.theta^2);
            % CovariogramModelChoice 3: parameters=[theta,p,sigma,sigmaError] covariance = sigmaError + obj.sigma^2*exp(-1/2*sum(distance^obj.p),obj.theta.^2)));
            % CovariogramModelChoice 4: parameters=[theta,sigma,sigmaError] covariance = sigmaError + obj.sigma^2*exp(-1/2*sum(distance^2),obj.theta.^2)));
            % CovariogramModelChoice 5: parameters=[theta,sigma,sigmaError] covariance = sigmaError + obj.sigma^2*(1 + sqrt(3)*weightedEuclidean).*exp((-sqrt(3)*weightedEuclidean)), 
            %                           weightedEuclidean=sqrt(distance.^2/obj.theta^2);
            % CovariogramModelChoice 6: as  CovariogramModelChoice 5 but with
            %                           weightedEuclidean = sqrt(sum(distance^2/obj.theta.^2))));
            % CovariogramModelChoice 7: covariance = sigmaError + obj.sigma^2*(1 + sqrt(5)*weightedEuclidean+5/3*weightedEuclidean.^2).*exp((- sqrt(5)*weightedEuclidean));
            %                           weightedEuclidean=sqrt(distance.^2/obj.theta^2);
            % CovariogramModelChoice 8: as  CovariogramModelChoice 7 but with
            %                           weightedEuclidean = sqrt(sum(distance^2/obj.theta.^2))));
            
            % Choose between the differnt models for the variogram
            Model = varargin{1};
            obj.CovariogramModelChoice = Model;
            
            % Adjust the covariogram parameter bounds (lower/upper) and
            % the first values
            obj.setCovariogramParameterBoundsAndValues(Model);
            
            obj.estimateVariance
        end

        %------------------------------------------------------------------
        function []=setCovariogramParameterBoundsAndValues(obj,varargin) 
            % Choose between the differnt models for the variogram
            % varargin = 1 ... exponential (sigma*exp(-(x/theta)^p))
            % varargin = 2 ... exponential with nugget (default) (sigmaError + sigma*exp(-(x/theta)^p))
            % varargin = 3 ... exponential with separated thethas and ps (multigaussian) (sigmaError + sigma*exp(-sum_(h=1)^k(theta(h)*abs(x(i,h)-x(j,h))^p(h)))
            Model = varargin{1};
            minValue = 1e-10;
            switch Model;
                case 1
                    obj.LBCovariogramModelParameters = [minValue,minValue];
                    obj.UBCovariogramModelParameters = [1e10,1e10];
                    obj.InitialCovariogramParameters = [];
                    obj.theta = 1;
                    obj.sigma = 1;
                    obj.nCovariogramParameters = 2;
                case {2,5,7}
                    obj.LBCovariogramModelParameters = [minValue,minValue,minValue];
                    obj.UBCovariogramModelParameters = [1e10,1e10,1e10];
                    obj.InitialCovariogramParameters = [];
                    obj.theta = 1;
                    obj.sigma = 1;
                    obj.sigmaError= 1;
                    obj.nCovariogramParameters = 3;
                case 3
                    if obj.nInputVar==0
                        error('nInputVar=0, please set InputData first');
                    end
                    switch Model
                        case 3
                            nTotal = obj.nInputVar;
                        case 4
                            nTotal = obj.nInputVar+nchoosek(obj.nInputVar,2);
                        otherwise
                    end
                    if obj.NormOutput
                        obj.LBCovariogramModelParameters = [ones(1,nTotal)*minValue,ones(1,nTotal)*minValue,minValue,minValue];
                        obj.UBCovariogramModelParameters = [ones(1,nTotal)*inf,ones(1,nTotal)*2,1,1];
                        obj.InitialCovariogramParameters = [];
                    else
                        obj.LBCovariogramModelParameters = [ones(1,nTotal)*minValue,ones(1,nTotal)*minValue,minValue,minValue];
                        obj.UBCovariogramModelParameters = [ones(1,nTotal)*inf,ones(1,nTotal)*2,inf,inf];
                        obj.InitialCovariogramParameters = [];
                    end
                    obj.theta = ones(1,nTotal);
                    obj.p = ones(1,nTotal);
                    obj.sigma = 1;
                    obj.sigmaError= 1;
                    obj.nCovariogramParameters = 2*nTotal+2;
                case {4,6,8}
                    if obj.nInputVar==0
                        error('nInputVar=0, please set InputData first');
                    end
                    
                    nTotal = obj.nInputVar;
                    if obj.NormOutput
                        obj.LBCovariogramModelParameters = [ones(1,nTotal)*minValue,minValue,minValue];
                        obj.UBCovariogramModelParameters = [ones(1,nTotal)*inf,1,1];
                        obj.InitialCovariogramParameters = [];
                    else
                        obj.LBCovariogramModelParameters = [ones(1,nTotal)*minValue,minValue,minValue];
                        obj.UBCovariogramModelParameters = [ones(1,nTotal)*inf,inf,inf];
                        obj.InitialCovariogramParameters = [];
                    end
                    obj.theta = ones(1,nTotal);
                    obj.sigma = 1;
                    obj.sigmaError= 1;
                    obj.nCovariogramParameters = nTotal+2;
                otherwise
                    error('Non-acceptable model choice. The parameter CovariogramModelChoice = %i is not allowed',Model);
            end
            
            % Decide for the distance calculation
            switch obj.CovariogramModelChoice
                case {1,2,5,7}
                    obj.CovariogramUsesEuclideanDistance = true;
                    obj.CovariogramUsesAbsoluteDistance = false;
                case {3,4,6,8}
                    obj.CovariogramUsesAbsoluteDistance = true;
                    obj.CovariogramUsesEuclideanDistance = false;
                otherwise
                    error('CovariogramModelChoice = %d is not defined',obj.CovariogramModelChoice)
            end
        end
        %------------------------------------------------------------------
        function []=setCR(obj,CR)
            % Set the fraction of the population at the next generation, not including elite children, that is created by the crossover function
            obj.CR = CR;
        end
        %------------------------------------------------------------------
        function []=setGenerations(obj,Generations)
            if Generations<0
                error('Generations must be positive')
            end
            % Set the Positive integer specifying the maximum number of iterations before the algorithm halts
            obj.Generations = Generations;
        end
        %------------------------------------------------------------------
        function []=setnIterationsSolver(obj,nIterationsSolver)
            if nIterationsSolver<0
                error('Iterations must be positive')
            end
            % Set the Positive integer specifying the maximum number of iterations before the algorithm halts
            obj.nIterationsSolver = nIterationsSolver;
        end
%         %------------------------------------------------------------------
%         function []=setMigrationFraction(obj,MigrationFraction)
%             % Set the Scalar between 0 and 1 specifying the fraction of
%             % individuals in each subpopulation that migrates to a
%             % different subpopulation
%             obj.MigrationFraction = MigrationFraction;
%         end
        %------------------------------------------------------------------
        function []=setNormInput(obj,NormInput)
            % Decide if the input shall be normalized(1) or not (0)
            if ~isempty(obj.InputData)
                warning('"InputData" is not empty. Please define "InputData" only after setting the "NormInput"-option')
            end
            obj.NormInput = NormInput;
        end
        %------------------------------------------------------------------
        function []=setNormOutput(obj,NormOutput)
            % Decide if the output shall be normalized(1) or not (0)
            if ~isempty(obj.OutputData)
                warning('"OutputData" is not empty. Please define "OutputData" only after setting the "NormOutput"-option')
            end
            obj.NormOutput = NormOutput;
%             error('Normalization of Output ist not allowed')
        end
        %------------------------------------------------------------------
        function []=setcheckVariogram(obj,checkVariogram)
            % (1) if the Output is normalized. (0) if not
            obj.checkVariogram = checkVariogram;
        end
        %------------------------------------------------------------------
        function []=setUseInverse(obj,UseInverse)
            % Decide if the inverse of the covariance matrix should be
            % calculated explicit (UseInverse=1)
            obj.UseInverse = UseInverse;
        end
        %------------------------------------------------------------------
        function []=setTimeLimit(obj,TimeLimit)
            if TimeLimit<0
                error('TimeLimit must be positive')
            end
            % Define how long the search algorithm shall maximal run
            obj.TimeLimit = TimeLimit ;
        end
        %------------------------------------------------------------------
        function []=setMeanCoeffBasisFct(obj,MeanCoeffBasisFct)
            % Define mean values of the coefficients of the basis functions used in universal kriging and bayesian kriging
            
            if  size(MeanCoeffBasisFct,1)~=obj.nBasisFct&&size(MeanCoeffBasisFct,2)~=obj.nBasisFct
                error('MeanCoeffBasisFct should be %ix1 vector but it is a %ix%i vector.',obj.nBasisFct,size(MeanCoeffBasisFct,1),size(MeanCoeffBasisFct,2));
            end
            if  size(MeanCoeffBasisFct,1)~=obj.nBasisFct&&size(MeanCoeffBasisFct,1)~=1
                error('MeanCoeffBasisFct should be %ix1 vector but it is a %ix%i vector.',obj.nBasisFct,size(MeanCoeffBasisFct,1),size(MeanCoeffBasisFct,2));
            end
            if  size(MeanCoeffBasisFct,1)~=obj.nBasisFct&&size(MeanCoeffBasisFct,2)==obj.nBasisFct
                MeanCoeffBasisFct = MeanCoeffBasisFct';
                warning('MeanCoeffBasisFct should be %ix1 vector but it is a %ix%i vector. THe vector is transposed',obj.nBasisFct,size(MeanCoeffBasisFct,1),size(MeanCoeffBasisFct,2));
            end 
            obj.MeanCoeffBasisFct = MeanCoeffBasisFct;
            
        end
        %------------------------------------------------------------------
        function []=setCoVarCoeffBasisFct(obj,CoVarCoeffBasisFct)
            % Define covariance matrix of the coefficients of the basis functions used in universal kriging and bayesian kriging
            obj.CoVarCoeffBasisFct = CoVarCoeffBasisFct;
        end       
        %------------------------------------------------------------------
        function []=setAverageProductOfOutputs(obj,AverageProductOfOutputs)
            % Define covariance matrix of the coefficients of the basis functions used in universal kriging and bayesian kriging
            obj.AverageProductOfOutputs = AverageProductOfOutputs;
        end
        %------------------------------------------------------------------
        function []=setProdWeightedFunctions(obj,ProdWeightedFunctions)
            % Define covariance matrix of the coefficients of the basis functions used in universal kriging and bayesian kriging
            obj.ProdWeightedFunctions = ProdWeightedFunctions;
        end
        %------------------------------------------------------------------
        function []=setSumSumOfWeightedBasisFcts(obj,SumSumOfWeightedBasisFcts)
            % Define covariance matrix of the coefficients of the basis functions used in universal kriging and bayesian kriging
            obj.SumSumOfWeightedBasisFcts = SumSumOfWeightedBasisFcts;
        end
        %------------------------------------------------------------------
        function []=setUseSolver(obj,UseSolver)
            obj.UseSolver = UseSolver;
        end
        %------------------------------------------------------------------
        function []=setCorrelationMatrix(obj,CorrelationMatrix)
            obj.CorrelationMatrix = CorrelationMatrix;
        end
        %------------------------------------------------------------------
        function []=setnBasisFctParameters(obj,nBasisFctParameters)
            obj.nBasisFctParameters = nBasisFctParameters;
        end
        %------------------------------------------------------------------
        function []=setBasisFctType(obj,BasisFctType)
            obj.BasisFctType = BasisFctType;
        end
        %------------------------------------------------------------------
        function []=setInitialBasisFctParameters(obj,InitialBasisFctParameters)
            if (length(InitialBasisFctParameters)>obj.nBasisFctParameters)
                warning('The number of initial parameter values(%d) does not match the number of parameters in the defined basis function(%d).\nMaybe you did not define the basis function yet. It is recommented to define them first and repeat the "setInitialBasisFctParameters" afterwards',length(InitialBasisFctParameters),obj.nBasisFctParameters)
            end
            if (length(InitialBasisFctParameters)<obj.nBasisFctParameters)
                error('The number of initial parameter values(%d) does not match the number of parameters in the defined basis function(%d).\nMaybe you did not define the basis function yet. It is recommented to define them first and repeat the "setInitialBasisFctParameters" afterwards',length(InitialBasisFctParameters),obj.nBasisFctParameters)
            end
            
            if (size(InitialBasisFctParameters,2)~=1&&size(InitialBasisFctParameters,1)~=1)
                error('Input should be a column vector of size %d',obj.nBasisFctParameters);
            elseif size(InitialBasisFctParameters,2)~=1&&size(InitialBasisFctParameters,1)==1
                warning('Input should be a column vector of size %d and not a row vector. The vector is transposed',obj.nBasisFctParameters)
                InitialBasisFctParameters = InitialBasisFctParameters';
            end
            obj.InitialBasisFctParameters = InitialBasisFctParameters;
        end
        %------------------------------------------------------------------
        function []=setBasisFctParameters(obj,BasisFctParameters)
            if (length(BasisFctParameters)>obj.nBasisFctParameters)
                warning('The number of initial parameter values(%d) does not match the number of parameters in the defined basis function(%d).\n',length(BasisFctParameters),obj.nBasisFctParameters)
            end
            
            if (length(BasisFctParameters)<obj.nBasisFctParameters)
                error('The number of initial parameter values(%d) does not match the number of parameters in the defined basis function(%d).\n',length(BasisFctParameters),obj.nBasisFctParameters)
            end
            obj.BasisFctParameters= BasisFctParameters;
            
            if ~isempty(obj.getCovariogramMatrix)
                obj.calcCovariogramMatrix
            end
        end
        %------------------------------------------------------------------
        function []=setLBBasisFctParameters(obj,LBBasisFctParameters)
            if length(LBBasisFctParameters)~=obj.nBasisFctParameters&&~isempty(LBBasisFctParameters)
                error('"LBBasisFctParameters" has to be nBasisParameter x 1')
            end
            obj.LBBasisFctParameters = LBBasisFctParameters;
        end
        %------------------------------------------------------------------
        function []=setUBBasisFctParameters(obj,UBBasisFctParameters)
            if length(UBBasisFctParameters)~=obj.nBasisFctParameters&&~isempty(UBBasisFctParameters)
                error('"UBBasisFctParameters" has to be nBasisParameter x 1')
            end
            obj.UBBasisFctParameters = UBBasisFctParameters;
        end
        %------------------------------------------------------------------
        function []=setPopulationSize(obj,PopulationSize)
            if PopulationSize<0
                error('PopulationSize must be positive')
            end
            obj.PopulationSize = PopulationSize;
        end
        %------------------------------------------------------------------
        function []=setShowDetails(obj,ShowDetails)
            if ~islogical(ShowDetails)
                error('ShowDetails should be boolean')
            end
            obj.ShowDetails = ShowDetails;
        end
        %------------------------------------------------------------------
        function []=setnIterationCrossValidation(obj,nIterationCrossValidation)
            obj.nIterationCrossValidation = nIterationCrossValidation;
        end
        %------------------------------------------------------------------
        function []=setnLowerUpperBoundSections(obj,nLowerUpperBoundSections)
            obj.nLowerUpperBoundSections = nLowerUpperBoundSections;
        end
        %------------------------------------------------------------------
        function []=setnTimeSections(obj,nTimeSections)
            obj.nTimeSections = nTimeSections;
        end
        %------------------------------------------------------------------
        function []=setnCrossOverSections(obj,nCrossOverSections)
            obj.nCrossOverSections = nCrossOverSections;
        end
        %------------------------------------------------------------------
        function []=setnCrossoverFractionSections(obj,nCrossoverFractionSections)
            obj.nCrossoverFractionSections = nCrossoverFractionSections;
        end
        %------------------------------------------------------------------
        function []=setnGenerationSections(obj,nGenerationSections)
            obj.nGenerationSections = nGenerationSections;
        end
        %------------------------------------------------------------------
        function []=setCrossValidationType(obj,CrossValidationType)
            obj.CrossValidationType = CrossValidationType;
        end
        %------------------------------------------------------------------
%         function []=setPlotCrossValidationResults(obj,PlotCrossValidationResults)
%             obj.PlotCrossValidationResults = PlotCrossValidationResults;
%         end
        %------------------------------------------------------------------
        function []=setCovariogramEstimationType(obj,CovariogramEstimationType)
            obj.CovariogramEstimationType= CovariogramEstimationType;
        end
        %------------------------------------------------------------------
        function []=setBasisFctCoefficients(obj,BasisFctCoefficients)
            obj.BasisFctCoefficients = BasisFctCoefficients;
        end
        %------------------------------------------------------------------
        function []=setInitialCovariogramParameters(obj,varargin)
            % [] = setUBCovariogramModelParameters(varargin)
            % Notice this function shall only be called AFTER
            % CovariogramModelChoice 1: varargin=[theta,sigma] (sigma*exp(-(x/theta)^p))
            % CovariogramModelChoice 2: varargin=[theta,sigma,sigmaError] (sigmaError + sigma*exp(-(x/theta)^p))
            % CovariogramModelChoice 3: varargin=[theta,sigma,p,sigmaError] (sigmaError + thetha*exp(-sum((x/sigma(i))^p(i)))
            % CovariogramModelChoice 3´4: varargin=[theta,sigma,sigmaError] (sigmaError + thetha*exp(-sum((x/sigma(i))^2))
            
            obj.checkIfNumberOfParameterValuesIsCorrect(varargin{1})
            
            
            if size(varargin{1},1)>1&&size(varargin{1},2)==1
                warning('InitialCovariogramParameters input should be a column vector. InitialCovariogramParameters is transposed')
                varargin{1} = varargin{1}';
            end
            
            obj.InitialCovariogramParameters = varargin{1};
        end
        %------------------------------------------------------------------
        function []=setOrdinaryRegressionKrigingObject(obj,OrdinaryRegressionKrigingObject)
            obj.OrdinaryRegressionKrigingObject = OrdinaryRegressionKrigingObject;
        end
        %------------------------------------------------------------------
        function []=setEstimatedOutputDataViaOrdinaryKriging(obj,EstimatedOutputDataViaOrdinaryKriging)
            obj.EstimatedOutputDataViaOrdinaryKriging = EstimatedOutputDataViaOrdinaryKriging;
        end
        %------------------------------------------------------------------
        function []=setAbsoluteInterpolation(obj,AbsoluteInterpolation)
            obj.AbsoluteInterpolation = AbsoluteInterpolation;
        end
        %------------------------------------------------------------------
        function []=setMaxSizeOfPredictions(obj,MaxSizeOfPredictions)
            obj.MaxSizeOfPredictions = MaxSizeOfPredictions;
        end
        %------------------------------------------------------------------
        function []=setRelativePartOfCrossValidation(obj,RelativePartOfCrossValidation)
            obj.RelativePartOfCrossValidation= RelativePartOfCrossValidation;
        end
        %------------------------------------------------------------------
        function []=setRelativePartOfConfidenceInterval(obj,RelativePartOfConfidenceInterval)
            obj.RelativePartOfConfidenceInterval = RelativePartOfConfidenceInterval;
        end
        %------------------------------------------------------------------
        function []=setShowWaitingBar(obj,ShowWaitingBar)
            if ~islogical(ShowWaitingBar)
                error('ShowWaitingBar should be boolean')
            end
            obj.ShowWaitingBar = ShowWaitingBar;
        end
        %------------------------------------------------------------------
        function []=setUseUniformDistributedPopulation(obj,UseUniformDistributedPopulation)
            obj.UseUniformDistributedPopulation = UseUniformDistributedPopulation;
        end
        %------------------------------------------------------------------
        function []=setCoeffLSQVargramCrossValCumLSQ(obj,CoeffLSQVargramCrossValCumLSQ)
            if(size(CoeffLSQVargramCrossValCumLSQ,1)==1&&size(CoeffLSQVargramCrossValCumLSQ,2)==3)
                obj.CoeffLSQVargramCrossValCumLSQ= CoeffLSQVargramCrossValCumLSQ;
            elseif(size(CoeffLSQVargramCrossValCumLSQ,1)==3&&size(CoeffLSQVargramCrossValCumLSQ,2)==1)
            obj.CoeffLSQVargramCrossValCumLSQ= CoeffLSQVargramCrossValCumLSQ';
            else 
                error('SIze of CoeffLSQVargramCrossValCumLSQ should 1x3 but it is %ix%i',size(CoeffLSQVargramCrossValCumLSQ,1),size(CoeffLSQVargramCrossValCumLSQ,2))
            end
        end
        %------------------------------------------------------------------
        function []=setPopInitMatrix(obj,PopInitMatrix)
            if(size(PopInitMatrix,2)~=2)
                error('PopInitMatrix must of size nPara x 2 but it is %ix%i',size(PopInitMatrix,1),size(PopInitMatrix,2))
            end
            if(sum(PopInitMatrix(:,1)>PopInitMatrix(:,2))>0)
                error('First column in PopInitMatrix has to contain the minimum values')
            end
            
            obj.PopInitMatrix = PopInitMatrix;
        end
        %------------------------------------------------------------------
        function []=setUseSimpleKriging(obj,UseSimpleKriging)
            obj.UseSimpleKriging = UseSimpleKriging;
        end
        %------------------------------------------------------------------
        function []=setUseParallel(obj,UseParallel)
            obj.UseParallel = UseParallel;
        end
        %------------------------------------------------------------------
        function []=setnCovariogramParameters(obj,nCovariogramParameters)
            if sum(size(nCovariogramParameters)>1)>0
                error('nCovariogramParameters shoudl be a scalar')
            end
            obj.nCovariogramParameters= nCovariogramParameters;
        end
        %------------------------------------------------------------------
        function []=setUseMatlabRegressionGP(obj,UseMatlabRegressionGP)
            necessaryFeatures = {'Statistics_Toolbox'};
            validLincence = cellfun(@(f) license('checkout',f),necessaryFeatures);
            versionStr = version;
            
            if ~islogical(UseMatlabRegressionGP)
                error('Inputput should be logical')
            end
            
            if UseMatlabRegressionGP&&(~validLincence||(str2double(versionStr(end-5:end-2))<2015)||(str2double(versionStr(end-5:end-2))==2015&&strcmp(versionStr(end-1),'a')))
                error('For Using RegressionGP, you need at least Matlab version 2015b and the Statistic Toolbox')
            end
            
            
            obj.UseMatlabRegressionGP = UseMatlabRegressionGP;
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
