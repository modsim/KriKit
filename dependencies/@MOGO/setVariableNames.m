function [] = setVariableNames(obj,inputNames,OutputNames)

    
    nObj = obj.KrigingAnalyzeObj.getnKrigingObjects;
    for iKrigObj = 1:nObj
        obj.KrigingAnalyzeObj.setInputVarNames(iKrigObj,inputNames);
        obj.KrigingAnalyzeObj.setKrigingObjectNames(iKrigObj,OutputNames{iKrigObj})
    end

end

