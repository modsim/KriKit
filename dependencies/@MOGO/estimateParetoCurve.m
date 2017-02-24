function [] = estimateParetoCurve(obj,iKrigingEstimation)
    nObj = obj.KrigingAnalyzeObj.getnKrigingObjects;

    obj.KrigingAnalyzeObj.predictParetoCurve(1:nObj,obj.nRealizations,(10^3),30);
    globalUncertainity=obj.KrigingAnalyzeObj.getGlobalParetoUncertainity;
    globalUncertainityNorm=obj.KrigingAnalyzeObj.getGlobalParetoUncertainityNorm;
    
    refPoint = (1-obj.KrigingAnalyzeObj.getMinMax(1:nObj))/2;
    transformedOutput = -bsxfun(@times,obj.KrigingAnalyzeObj.getMinMax(1:nObj),obj.OutputData);
    HV=Hypervolume_MEX(transformedOutput,refPoint);
    
    obj.GlobalUncertainity(iKrigingEstimation) = globalUncertainity;
    obj.GlobalUncertainityNorm(iKrigingEstimation) = globalUncertainityNorm;
    obj.HVVec(iKrigingEstimation) = HV;
    
    
end

