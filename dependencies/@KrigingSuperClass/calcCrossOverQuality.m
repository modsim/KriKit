function [quality]=calcCrossOverQuality(obj,varargin)
%[quality]=calcCrossOverQuality(obj,covariogramModelParameters)
%
% This function calculated the quality of test values of the coavriogram
% parameters by cross validation. I.e., for each sample location, the data
% point is removed and the Kriging prediction is calculated. Beased on
% both, the Kriging estimation of the expected value and the standard
% deviation, the quality of the given covariogram paraset is determined.
% The qualityi s defined based on the different Cross-validation types:
% 1 ... (1-1/nExp*sum(((output_true-output_predicted)/sigma_predicted)^2))^2
% 2 ... sum((output_true-output_predicted)/^2)
% 3 ... sum(((output_true-output_predicted)/sigma_predicted)^2
% 4 ... Combines option 2 and 3 by  
%       qualityTotal = quality1 + quality2/nUniqueInput/max(OutputBackup)^2;
%       modification of quality2 makes sure that both qualtiy have the same
%       scaling range/order of magnitude
% By default 4
%
% In the case of multiple measurements all data at the same point are
% erased and replaced by the mean values at this point
%
% Input: 
% covariogramModelParameters ... 1XnCovariogramParameters, vector which
%                                contains the test values of the covariance
%                                parameters
%
% You can set:
% - CrossValidationType ... decide which cross validation qulity definition
%                           is used
% - CovariogramModelChoice ... defined which covariogram is used
%
% You can get:
%
% Note: The following member variables are temporary changed. If you cancel
% the process they might have not returned to thier original value:
% - CovariogramModelParameters;
% - InputData;
% - OutputData;
% - InputData;
% - OutputData;
% - nExperiments;  
% - Inverse;
% - CovariogramMatrix;
%
% For imporivng performance. The original inverse of the covariance matrix
% is only updated, for each data poin which is removed. Algorithm for
% updating the inverse was implemented by (HAGER, 1989,'UPDATING THE
% INVERSE OF A MATRIX)
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

    % Initialization: Backup everything that is temporary change in this
    % function
    covParaBackup = obj.getCovariogramModelParameters;
    InputBackup = obj.getInputData;
    OutputBackup = obj.getOutputData;
    InputBackupNorm = obj.InputData;
    OutputBackupNorm = obj.OutputData;
    nExpBackup = obj.nExperiments;  
    UseInverseBackUP = obj.UseInverse;
        % Set the defined model parameters
    obj.setCovariogramModelParameters(varargin{1});
    covariogramBackupMatrix = obj.getCovariogramMatrix;
    
    % Calculated initial inverse of the covariogramgram matrix. Later this
    % inverse is only updated
    obj.UseInverse = true;
    inverseBackUP = inv(covariogramBackupMatrix);

    % Calculate the quality
    quality = 0;
    quality2 = 0;
    
    % In the case of multiple measurements all data at the same point are
    % erased and replaced by the mean values at this point
    uniqueInput = unique(obj.getInputData,'rows');
    
    
    for iExpOut = 1:size(uniqueInput,1)
        
        % Reset input and output data
        InputData = InputBackup;
        InputDataNorm = InputBackupNorm;
        OutputData = OutputBackup;
        OutputDataNorm = OutputBackupNorm;
        
        % Erase all rows and column which are associated with current
        % data point
            % Find all Experiments using the same input parameters
        [~,binVector]=ismember(InputBackup,uniqueInput(iExpOut,:),'rows');
        indexSameInput = find(binVector==1);
            % Erase
        InputData(indexSameInput,:)         = [];
        InputDataNorm(indexSameInput,:)     = [];
        OutputData(indexSameInput,:)        = [];
        OutputDataNorm(indexSameInput,:)    = [];
            % Replaceoriginal data set
        setInputData(InputData,InputDataNorm);
        setOutputData(OutputData,OutputDataNorm);
        
        % Update covariance matrix (erase lines which are associated with
        % erased data set)
        indicesKeep = setdiff(1:size(covariogramBackupMatrix,1),indexSameInput);
        covariogramMatrix = covariogramBackupMatrix(indicesKeep,indicesKeep);
        obj.CovariogramMatrix = covariogramMatrix;
        
        % Instead of calculating the inverse each time de novo, the old
        % inverse is just updated (HAGER, 1989,'UPDATING THE INVERSE OF A
        % MATRIX*'). Only possible in step for when no multiple
        % measurements exists at the point of current investigation.
        if length(indexSameInput)==1
            obj.InvCovariogramMatrix=invupdatered(inverseBackUP,indexSameInput,indexSameInput);
        else
            obj.InvCovariogramMatrix = inv(covariogramMatrix);
        end

        % Make the Kriging prediction
        [output]=obj.prediction(uniqueInput(iExpOut,:));
        expected = output(:,1);
        sigma = output(:,2);
        
        % Compare wiht the true value
        switch obj.CrossValidationType
            case {1,3}
                proto = ((mean(OutputBackup(indexSameInput))-expected)/sigma)^2;

                if isnan(proto)
                    quality = realmax;
                else
                    quality = quality + proto;
                end
            case 2
                if isnan(expected)
                    quality = realmax;
                else
                    quality = quality + (mean(OutputBackup(indexSameInput))-expected)^2;
                end
            case 4
                proto = ((mean(OutputBackup(indexSameInput))-expected)/sigma)^2;

                if isnan(proto)
                    quality = realmax;
                else
                    quality = quality + proto;
                end
                    quality2 = quality2 + (mean(OutputBackup(indexSameInput))-expected)^2;
            otherwise
                error('CrossValidationType: %d is not defined\n',obj.CrossValidationType)
        end
    end
    
    switch obj.CrossValidationType
        case 1
            quality = abs(1-quality/size(uniqueInput,1));
        case 4
            quality = abs(1-quality/size(uniqueInput,1));
%             quality = quality*quality2;
            quality = quality+quality2/size(uniqueInput,1)/max(OutputBackup.^2);
        otherwise
                
    end
    
    % Restore old values
    obj.checkVariogram = 0;
    obj.setInputData(InputBackup);
    obj.setOutputData(OutputBackup);
    obj.CovariogramMatrix = covariogramBackupMatrix;
    obj.InvCovariogramMatrix = inverseBackUP;
    obj.nExperiments = nExpBackup;
    obj.UseInverse = UseInverseBackUP;
    obj.setCovariogramModelParameters(covParaBackup);
    
    %% Set Functions
    % ---------------------------------------------------------------------
    function [] = setInputData(Input,InputNorm)
        % Replace current input data set by modify input
        
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
                obj.InputData = InputNorm;
            otherwise
                error('NormInput = %i is not defined or Crossvalidation',obj.NormInput);
        end
    end
    % ---------------------------------------------------------------------
    function [] = setOutputData(Output,OutputNorm)
        % Replace current output data set by modify output
        
        obj.OutputData = Output;
        obj.OutputData_True = Output;
        obj.nOutputVar = obj.nExperiments;

        switch obj.NormOutput
            case 0

            case 1
                obj.OutputData = OutputNorm;
            otherwise
                error('NormOutput = %i is not defined for cross validation',obj.NormOutput);
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
