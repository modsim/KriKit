classdef MOGO<handle
    
    %% Public Members
    properties(GetAccess='public',SetAccess='public')
       KrigingAnalyzeObj = [];
       CovModelChoice = 6;
       SaveInputDataOverIterations = cell(0);
       SaveOutputDataOverIterations = cell(0);
       SaveCovarParametersOverIterations = cell(0);
       InputDataRange = [];
       nInputVar = 0;
       % Number of realizations drawn during conditional simulation
       nRealizations = 5e2;
       
       InputData = [];
       OutputData = [];
    end
    properties(GetAccess='public',SetAccess='protected')
       
       % Saving optimization results in each iteration
       GlobalUncertainity = [];
       GlobalUncertainityNorm = [];
       HVVec = [];
    end
    %% Private Members
    properties(GetAccess='private',SetAccess='private')

    end
    
    %% Protected Members
    properties(GetAccess='protected',SetAccess='protected')
        
    end
    
    methods(Abstract)
        % These function have to be user defined in each inherinted class
        
        % ----------------------------------------------------------------
        % Used in "generateData"
        [output] = objFct(obj,varargin)
    end
    
    methods
        %% Constructor
        function obj = MOGO()
            obj.KrigingAnalyzeObj = BayesianOptimizationClass();
            obj.KrigingAnalyzeObj.addKrigingObject(1,'Test Name')
            obj.KrigingAnalyzeObj.KrigingObjects{1}.setNormInput(true)
            obj.KrigingAnalyzeObj.KrigingObjects{1}.setNormOutput(true)
            obj.KrigingAnalyzeObj.KrigingObjects{1}.setShowDetails(true)
        end
        %% Copy Operator for a shallow copy
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
        %% General Methods
        [] = calcKriging(obj,varargin)
        % ----------------------------------------------------------------
        function []=saveObject(obj,varargin)
            index = varargin{1};
            save(strcat('tmpSave_',num2str(index)),'obj')
        end
        % ----------------------------------------------------------------
        [] = saveTmpResults(obj,iIter)
        % ----------------------------------------------------------------
        [] = estimateParetoCurve(obj,iKrigingEstimation)
        % ----------------------------------------------------------------
        [] = setVariableNames(obj,inputNames,OutputNames)
        % ----------------------------------------------------------------
        [] = calcAndSaveNewSamples(obj,iterationNumber)
        % ----------------------------------------------------------------
        [inputData,outputData]=collectAllDataUntilGivenIteration(obj,finalIteration)
        % ----------------------------------------------------------------
        function []=SaveCurrentCovarParameter(obj,iterationNumber)
            obj.SaveCovarParametersOverIterations{iterationNumber} = obj.KrigingAnalyzeObj.KrigingObjects{1}.getCovariogramModelParameters;
        end
        % ----------------------------------------------------------------
        []=generateNewData(obj,varargin);
        %% Set function
        function []=setInputDataRange(obj,InputDataRange)
            % Check 
            if size(InputDataRange,1)== 2 && size(InputDataRange,2)== obj.nInputVar && obj.nInputVar~=2
                InputDataRange = InputDataRange';
            end
            obj.InputDataRange = InputDataRange;
        end
        
        function []=setInputData(obj,InputData)
            obj.resetData
            obj.InputData = InputData;
            obj.KrigingAnalyzeObj.KrigingObjects{1}.setInputData(obj.InputData)
            obj.SaveInputDataOverIterations{1} = obj.InputData;
            obj.nInputVar = size(InputData,2);
        end
        % ----------------------------------------------------------------
        function []=setOutputData(obj,OutputData)
            obj.OutputData = OutputData;
            obj.KrigingAnalyzeObj.KrigingObjects{1}.setOutputData(obj.OutputData)
            obj.SaveOutputDataOverIterations{1} = obj.OutputData;
        end
    end
end

