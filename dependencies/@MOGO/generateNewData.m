function [] = generateNewData(obj,varargin)
    iterationNumber = varargin{1};
    
    outputProto = obj.objFct(obj.SaveInputDataOverIterations{iterationNumber});
    obj.SaveOutputDataOverIterations{iterationNumber} = outputProto(:,1:end-obj.nInputVar);
        
    % Save Input Data
    obj.SaveInputDataOverIterations{iterationNumber} = outputProto(:,end-obj.nInputVar+1:end);
    
end

