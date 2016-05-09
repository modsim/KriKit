function [covariance] = CovarModel(obj,distance,varargin)
% [covariance] = CovarModel(distance,useSigmaError)
% This function use the chosen mode to calculate the covariance between two
% points as function of their distance. The distance can be a Matrix
% nSamplePoints x nInputVar.
% If useSigmaError=true. the covariogram model parameter sigmaError is
% considered for simulation of measurment error
%
%
% You can set: 
% - CovariogramModelChoice .. define which covariogram is used
% - Covariogram parameters ... theta, p, sigma, sigmaError
%
% You can get: -
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.


% Note For Programmer: The following lines make sure that sigmaError only
% plays a role when ignoreSigmaError=true
% Euclidean Distance:
% if ~isempty(varargin)&&varargin{1}==1
%     [r,c] = find(distance==0);
%     covariance(r,c) = covariance(r,c) + obj.sigmaError^2;
% end
%---- Or ---
% Absolute Distance:
% [covariance]=decideForOffset(covariance);
%
% If absolute Distance is calculate:
% Distance is saved in a tensor of size: 
% (nInputDataXnInputData)X(nInputVar,1)X(nInputDataXnInputData) and is interpreted as
% [pairwise distance of with respect to input variable 1]X[scalar for input variable index to which is comapred with]X[pairwise distance of with respect to input variable 2]

% Decided between differen covariogram models
switch obj.CovariogramModelChoice
    case 1
        covariance = obj.sigma^2*exp(-1/2*distance.^2/obj.theta^2);
    case 2
        
        covariance = obj.sigma^2*exp(-1/2*distance.^2/obj.theta^2);
        
        if ~isempty(varargin)&&varargin{1}==1
            [r,c] = find(distance==0);
            covariance(r,c) = covariance(r,c) + obj.sigmaError^2;
        end
    case {3,4}
        % Calculate the distances
        switch obj.CovariogramModelChoice
            case 3
                covariance = obj.sigma^2*exp(-1/2*sum(bsxfun(@rdivide,bsxfun(@power,distance,obj.p),obj.theta.^2),2));
            case 4
    %                     covarianceAs = obj.sigma^2*exp(-1/2*distance.^2/obj.theta^2);
                covariance = obj.sigma^2*exp(-1/2*sum(bsxfun(@rdivide,bsxfun(@power,distance,2),obj.theta.^2),2));
        end

        [covariance]=decideForOffset(covariance);
     case {5,6}
         switch obj.CovariogramModelChoice
             case 5
                weightedEuclidean = sqrt(distance.^2/obj.theta^2);
             case 6
                weightedEuclidean = sqrt(sum(bsxfun(@rdivide,bsxfun(@power,distance,2),obj.theta.^2),2));
         end

         covariance = obj.sigma^2*(1 + sqrt(3)*weightedEuclidean).*exp((- sqrt(3)*weightedEuclidean));

         switch obj.CovariogramModelChoice
             case 5
                if ~isempty(varargin)&&varargin{1}==1
                    [r,c] = find(distance==0);
                    covariance(r,c) = covariance(r,c) + obj.sigmaError^2;
                 end
             case 6
                [covariance]=decideForOffset(covariance);
         end
     case {7,8}
         switch obj.CovariogramModelChoice
             case 7
                weightedEuclidean = sqrt(distance.^2/obj.theta^2);
             case 8
                weightedEuclidean = sqrt(sum(bsxfun(@rdivide,bsxfun(@power,distance,2),obj.theta.^2),2));
         end

         covariance = obj.sigma^2*(1 + sqrt(5)*weightedEuclidean+5/3*weightedEuclidean.^2).*exp((- sqrt(5)*weightedEuclidean));

         switch obj.CovariogramModelChoice
             case 7
                if ~isempty(varargin)&&varargin{1}==1
                    [r,c] = find(distance==0);
                    covariance(r,c) = covariance(r,c) + obj.sigmaError^2;
                 end
             case 8
                [covariance]=decideForOffset(covariance);
         end
    otherwise
        error('Non-acceptable model choice. The parameter CovariogramModelChoice = %i is not allowed',obj.CovariogramModelChoice);
end
%% Nested Functions
function [covariance]=decideForOffset(covariance)
    if ~isempty(varargin)&&varargin{1}==1
        % It is easier when only two input variable exists
        if length(size(covariance))==2
            [iI]=find(ismember(distance,zeros(1,size(distance,2)),'rows'));
            covariance(iI) = covariance(iI) + obj.sigmaError^2;
        else
            covariance=calcCovarianceHigherDimensions(covariance);
        end
    end
end
% ------------------------------------------------------------------------
function [covariance]=calcCovarianceHigherDimensions(covariance)
    % Initialization
    iI_total = zeros(size(distance,3),1);
    iE_total = zeros(size(distance,1),1);
    nMultipleMeasurements = zeros(size(distance,3),1);

    indexI = 0;
    % Find all interpolation points which are identical to
    % the the provided measurement points (distance is for all input variable = 0)
    for iInputPoints = 1:size(distance,3)
        zeroRows = ones(size(distance,1),1);
        % Find indices of rows where all columns are equal to zero (faster than using find)
        for iInputVar = 1 : obj.nInputVar
            zeroRows = zeroRows.*(distance(:,iInputVar,iInputPoints)==0);
        end
        zeroRows = find(zeroRows==1);

        if ~isempty(zeroRows)
            % Save the rows which are zero and save the
            % input index
            iE_total(indexI+1:indexI+length(zeroRows)) = zeroRows;
            iI_total(indexI+1:indexI+length(zeroRows)) = iInputPoints;
            
            % In the case of multiple measurements the average
            % value is chosen at this point. Therefore each
            % measurement get the same amount of extra covariance
            nMultipleMeasurements(indexI+1:indexI+length(zeroRows)) = length(zeroRows);
            indexI=indexI+length(zeroRows);
        end
    end
    iE_total(iE_total<=0)=[];
    iI_total(iI_total<=0)=[];

    % Final Output
    for iChangeInCovarianceMatrix=1:length(iE_total)
        covariance(iE_total(iChangeInCovarianceMatrix),:,iI_total(iChangeInCovarianceMatrix)) = covariance(iE_total(iChangeInCovarianceMatrix),:,iI_total(iChangeInCovarianceMatrix)) + obj.sigmaError^2;
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
