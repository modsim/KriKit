function [Data,dataShown,OutlierShown]=determineRelevantDataPointsForPlot(obj,varargin)
% [Data,dataShown,OutlierShown]=determineRelevantDataPointsForPlot(obj,krigingObjIndex,dimensionInterpolation)
% 
%  This function determines which data points are relevant for the current
%  plot.
%
% Input:
% KrigingObjectIndex ... index of Kriging Object(s) used when
%                        "calcMutualInterpolation_23D()" was called
% dimensionInterpolation ... decision if 2D or 3D interpolation shall be
%                            created. dimensionInterpolation has to be
%                            consistent with the call of
%                            "calcMutualInterpolation_23D()". In case of
%                            screening dimensionInterpolation=4
%
% Output: 
% Data ... total set of sample points used for creating the kriging
%          model
% dataShown ... set of sample points which are relevant for the plot
% OutlierShown ... set of sample points which are considered to be outlier
%
%
% 
%  You can set: 
%  - ShowOutlier ... if true outlier are determined. For further details
%                    see "findPotentialOutlier()"
%
%  You can get: -
% 
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

iKrigingObject = varargin{1};
dimensionInterpolation = varargin{2};

Data = obj.KrigingObjects{iKrigingObject}.getInputData;
% Plot only data which actual met the desired values
    % Copy Data
Data2 = Data;
switch dimensionInterpolation
    case 2
        Data2(:,obj.KrigingPrediction_Interpolation2D{iKrigingObject,3})=[];

        % Don't consider the current investiagted variable
        InputData = obj.KrigingPrediction_Interpolation2D{iKrigingObject,2};
        InputData(:,obj.KrigingPrediction_Interpolation2D{iKrigingObject,3}) = [];
    case 3
        Data2(:,obj.KrigingPrediction_Interpolation3D{iKrigingObject,3})=[];

        % Don't consider the current investiagted variable
        InputData = obj.KrigingPrediction_Interpolation3D{iKrigingObject,2};
        InputData(:,obj.KrigingPrediction_Interpolation3D{iKrigingObject,3}) = [];
    otherwise
            error('Plotting function is only allowed for 2D and 3D interpolation')
end

% The remaining variable have only one value in each row
InputData = InputData(1,:); 

% Find Rows
% Consider only data in the defined region
isInInvestigatedSection = obj.findIndicesOfDataForPlot(iKrigingObject,Data);
r = ismember(Data2,InputData,'rows');
r = logical(r.*isInInvestigatedSection);


if obj.getShowOutlier
    obj.findPotentialOutlier(iKrigingObject);
    potentialOutlierBinary = obj.getPotentialOutlierBinary;
else
    potentialOutlierBinary = false(size(Data,1),1);
end
dataShown= r&(~potentialOutlierBinary);
OutlierShown = r&potentialOutlierBinary;

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
