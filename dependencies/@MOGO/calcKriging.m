function [] = calcKriging(obj,varargin)

    % Create kriging model for each object
    for iKrigingObject = 1:obj.KrigingAnalyzeObj.getnKrigingObjects
        obj.KrigingAnalyzeObj.KrigingObjects{iKrigingObject}.setUseMatlabRegressionGP(false);
        
        obj.KrigingAnalyzeObj.KrigingObjects{iKrigingObject}.setInputData(obj.InputData);
        obj.KrigingAnalyzeObj.KrigingObjects{iKrigingObject}.setOutputData(obj.OutputData(:,iKrigingObject));

        obj.KrigingAnalyzeObj.KrigingObjects{iKrigingObject}.setCovariogramModelChoice(obj.CovModelChoice)
        maxValue = 1e2;
        
        if obj.KrigingAnalyzeObj.KrigingObjects{iKrigingObject}.getCovariogramModelChoice==3
            obj.KrigingAnalyzeObj.KrigingObjects{iKrigingObject}.setLBCovariogramModelParameters([ones(1,obj.KrigingAnalyzeObj.KrigingObjects{1}.getnCovariogramParameters)]*1e-10)
            obj.KrigingAnalyzeObj.KrigingObjects{iKrigingObject}.setUBCovariogramModelParameters([maxValue,maxValue,2,2,maxValue,maxValue])
        else
            obj.KrigingAnalyzeObj.KrigingObjects{iKrigingObject}.setLBCovariogramModelParameters([ones(1,obj.KrigingAnalyzeObj.KrigingObjects{1}.getnCovariogramParameters)]*1e-10)
            obj.KrigingAnalyzeObj.KrigingObjects{iKrigingObject}.setUBCovariogramModelParameters([ones(1,obj.KrigingAnalyzeObj.KrigingObjects{1}.getnCovariogramParameters)]*maxValue)
        end
        
        
        obj.KrigingAnalyzeObj.KrigingObjects{iKrigingObject}.generateRegressionGPModel()
%         obj.KrigingAnalyzeObj.KrigingObjects{iKrigingObject}.calcCovariogramMatrix
        
    end
    
    obj.setInputDataRange(obj.InputDataRange(1:obj.nInputVar,:));
end