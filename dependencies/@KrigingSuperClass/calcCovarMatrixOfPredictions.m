function [] = calcCovarMatrixOfPredictions(obj,inputData)
% This function calculate the covariance matrix of the predictions. THis
% function should be run right after the Kriging prediction. Do not use
% this function if bayesian or disjunctive Kriging is used
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.



% Calculate all Combinations
nDataPoints = size(inputData,1);

if obj.NormInput==1
    for iInput = 1:obj.nInputVar
        inputData(:,iInput) = obj.scale(inputData(:,iInput), min(obj.InputData_True(:,iInput)),max(obj.InputData_True(:,iInput)),0,1);
    end
end
        

[nComb1,nComb2] = ndgrid(1:nDataPoints,1:nDataPoints);
nComb = [nComb1(:),nComb2(:)];
nOriginalDataPoints = obj.getnExperiments;
CovMatrixofExperiments = obj.getCovariogramMatrix;

% Check if all weights and covariance values are saved 
if size(obj.Weights,2)<max(nComb(:,1))||size(obj.CovargramVectors,2)<max(nComb(:,2))
    error('Not all weights are saved. Most likely the value of "MaxSizeOfPredictions" is low. It should be at least %i',max(max(nComb)))
end

% covariance_zero = ones(nDataPoints,nDataPoints)*obj.CovarModel(zeros(1,size(obj.DistInput,2)),1);
% obj.Weights = obj.Weights*(max(obj.getOutputData)-min(obj.getOutputData));

% Calculate the Covariance Matrix parts 1:nOriginalDataPoints
part1 = reshape(diag(obj.Weights(1:nOriginalDataPoints,nComb(:,1))'*obj.CovargramVectors(1:nOriginalDataPoints,nComb(:,2))),nDataPoints,nDataPoints);
part2 = reshape(diag(obj.Weights(1:nOriginalDataPoints,nComb(:,2))'*obj.CovargramVectors(1:nOriginalDataPoints,nComb(:,1))),nDataPoints,nDataPoints);
part3 = reshape(diag(obj.Weights(1:nOriginalDataPoints,nComb(:,1))'*CovMatrixofExperiments(1:nOriginalDataPoints,1:nOriginalDataPoints)*obj.Weights(1:nOriginalDataPoints,nComb(:,2)))...
                ,nDataPoints,nDataPoints);
% part4 = 


% All but the diagonals
if obj.getCovariogramUsesEuclideanDistance
    % Calculate the distances between measurements
    delta = sqrt(sum((bsxfun(@minus,inputData(nComb(:,1),:),inputData(nComb(:,2),:)).^2),2));
    delta = reshape(delta,sqrt(size(delta,1)),sqrt(size(delta,1)));
        % Calculate the variogram values for theses distances
    part4 =  obj.CovarModel(delta,0);
elseif obj.getCovariogramUsesAbsoluteDistance
    % Calculate the variogram values for the input distances
        delta=abs(bsxfun(@minus,inputData(nComb(:,1),:),inputData(nComb(:,2),:)));
        nRowsColumns=sqrt(size(delta,1));
        d = ones(nRowsColumns,nRowsColumns);

        for iN = 1:nRowsColumns
            for iV = 1:obj.nInputVar
               d(:,(iN-1)*obj.nInputVar+iV)=delta(nRowsColumns*(iN-1)+1:nRowsColumns*(iN),iV);
            end
        end
        
        delta = reshape(d,nRowsColumns,obj.nInputVar,nRowsColumns);
        part4 =  reshape(obj.CovarModel(delta,0),nRowsColumns,nRowsColumns);
else
    error('either CovariogramUsesEuclideanDistance or CovariogramUsesAbsoluteDistance has to be set')
end

% switch obj.CovariogramModelChoice
%     case {1,2,5,7}
%             % Calculate the distances between measurements
%         delta = sqrt(sum((bsxfun(@minus,inputData(nComb(:,1),:),inputData(nComb(:,2),:)).^2),2));
%         delta = reshape(delta,sqrt(size(delta,1)),sqrt(size(delta,1)));
%             % Calculate the variogram values for theses distances
%         part4 =  obj.CovarModel(delta,0);
%     case {3,4,6,8}
%         % Calculate the variogram values for the input distances
%         delta=abs(bsxfun(@minus,inputData(nComb(:,1),:),inputData(nComb(:,2),:)));
%         nRowsColumns=sqrt(size(delta,1));
%         d = ones(nRowsColumns,nRowsColumns);
% 
%         for iN = 1:nRowsColumns
%             for iV = 1:obj.nInputVar
%                d(:,(iN-1)*obj.nInputVar+iV)=delta(nRowsColumns*(iN-1)+1:nRowsColumns*(iN),iV);
%             end
%         end
%         
%         delta = reshape(d,nRowsColumns,obj.nInputVar,nRowsColumns);
%         part4 =  reshape(obj.CovarModel(delta,0),nRowsColumns,nRowsColumns);
%     otherwise
%             error('Invalide CovariogramModelChoice in calcCovarMatrixOfPredictions')
% end
    % Diagonal has to be 0
part4 = part4 + eye(size(part4,1))*obj.sigmaError^2;



% Final Result
obj.CovarMatrixOfPredictions = part4-part1-part2+part3;
if obj.getNormOutput
    obj.CovarMatrixOfPredictions = obj.CovarMatrixOfPredictions*(max(obj.getOutputData)-min(obj.getOutputData))^2;
end

% if obj.NormOutput
%     obj.CovarMatrixOfPredictions = (obj.CovarMatrixOfPredictions-0)*(max(obj.getOutputData)-min(obj.getOutputData))^2;
% end

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
