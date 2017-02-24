function [newSamplePoint] = calcAndSaveNewSamples(obj,varargin)
%     [] = calcAndSaveNewSamplesMultiObj(obj,iterationNumber,OptimizationAlgorithm,UseOnlyMaximumValue)
    iterationNumber = varargin{1};
    
    if obj.KrigingAnalyzeObj.getnKrigingObjects>1
        obj.KrigingAnalyzeObj.determineParetoSet(1:obj.KrigingAnalyzeObj.getnKrigingObjects)
    end
%     obj.nParetoSamples = obj.KrigingAnalyzeObj.getnParetoSetExperiments;
            
    % Save the covariogram parameters
    obj.SaveCovarParametersOverIterations{iterationNumber-1} = obj.KrigingAnalyzeObj.KrigingObjects{1}.getCovariogramModelParameters;
    newSamplePoint = obj.KrigingAnalyzeObj.calcNewSamplesViaMCMC(1:obj.KrigingAnalyzeObj.getnKrigingObjects,'DRAM');


    % Generate New data points (For ome cases input variables are aslo part of the ouput of the obj function)
    obj.SaveInputDataOverIterations{iterationNumber} = newSamplePoint;
    obj.generateNewData(iterationNumber);

    % Save data points
    obj.InputData = [obj.InputData;obj.SaveInputDataOverIterations{iterationNumber}];
    obj.OutputData = [obj.OutputData;obj.SaveOutputDataOverIterations{iterationNumber}];
%             obj.KrigingAnalyzeObj.getNewSamples
end