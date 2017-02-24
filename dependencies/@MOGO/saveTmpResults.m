function [] = saveTmpResults(obj,iIter)
    %% Backup
    fileNameUnc = strcat('globalUncertainityVec_iIter',num2str(iIter));
    fileNameIn = strcat('input_iIter',num2str(iIter));
    fileNameOut = strcat('output_iIter',num2str(iIter));
    fileNameNormUnc = strcat('globalUncertainityNormVec_iIter',num2str(iIter));
    fileNameHV = strcat('HVMatrixVec_iIter',num2str(iIter));    
    
    % Save
    input = obj.InputData;
    output = obj.OutputData;
    GU  = obj.GlobalUncertainity;
    GUN = obj.GlobalUncertainityNorm;
    HVVec = obj.HVVec;
    save(fileNameIn,'input')
    save(fileNameOut,'output')
    save(fileNameUnc,'GU')
    save(fileNameNormUnc,'GUN')
    save(fileNameHV,'HVVec')
    
end

