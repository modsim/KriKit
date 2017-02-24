function [inputData,outputData]=collectAllDataUntilGivenIteration(obj,finalIteration)
    % First allocation
    nData=0;
    for iIter=1:finalIteration
        nData = nData + size(obj.SaveInputDataOverIterations{iIter},1);
    end

    % Now collecting
    inputData = zeros(nData,obj.nInputVar);
    outputData = zeros(nData,1);
    iRow = 0;
    for iIter=1:finalIteration
        nDataIter = size(obj.SaveInputDataOverIterations{iIter},1);
        inputData(iRow+1:iRow+nDataIter,:) = obj.SaveInputDataOverIterations{iIter};
        outputData(iRow+1:iRow+nDataIter,1) = obj.SaveOutputDataOverIterations{iIter};
        iRow = iRow + nDataIter;
    end
end