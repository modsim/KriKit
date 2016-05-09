function []=calcAndPlotInterpolation_23D_BestChoice(obj,varargin)
% []=calcAndPlotInterpolation_23D_BestChoice(KrigingObjectIndex,inputIndices,dimensionInterpolation)
%
%
% Kriging interpolation in the 2D/3D space (two input variable + output
% variable). Value of the fixes in input variable are automatically
% adjusted. Best choice refers to setting the values of input which are
% not on the x/y-axis varied to most abundant values in the data set or
% input values where at least "nMinimumDataPoint" data point show up in the
% interpolation plot. 
% E.g.:
% InputData = [1,2;...
%              2,2;...
%              3,3;...
%              4,1];
% obj.setnMinimumDataPoint(-1);
% By calling calcAndPlotInterpolation_23D_BestChoice(1,1,2) a
% 2D-Interpolation Plot is created where first input variable is varied
% over the x-axis and the second input variable is hold at 2 this value is
% most abundant in the data set
%
% For further details about the Kriging interpolation see documentation of
% "calcMutualInterpolation_23D()" 
%
% Input: 
% KrigingObjectIndex ... Index of kriging objects of interest
% inputIndices ... Indices of input variable which shall be veid on the x-
%                  (and y-)axis
% dimensionInterpolation ... Decide if this is a 2D or 3D interpolation
%                            (dimensionInterpolation=2/3) 
%
% You can set:
% nMinimumDataPoint ... minimal number of data points which
%                       should be contained by the plot. If negative, the
%                       most abundant value in the data set is chosen.
%
% You can get: 
% ChosenCombinationsForPlot ... chosen set of input variable values
%
% Note: If settings lead to sevaral valid sets values for the non-varied
% input variables, several interpolation plots are generated
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

%% Initialization
KrigingObjectIndex = varargin{1};
KrigingObjectIndex = obj.checkKrigingIndizes(KrigingObjectIndex);
inputIndices = varargin{2};
dimensionInterpolation = varargin{3};

checkLengthOfIndices(dimensionInterpolation);

indices2 = setdiff(1:obj.KrigingObjects{KrigingObjectIndex(1)}.getnInputVar,inputIndices);

nData = 0;
for iKriging=1:length(KrigingObjectIndex)
    nData = nData+obj.KrigingObjects{iKriging}.getnExperiments;
end

%% Find appropriate combination of input variable choice
ParacombinationsMatrix = [];
ambundanceVector = [];
for iKriging=KrigingObjectIndex'
    [Paracombinations,ambundance] = obj.getUniqueInputVarCombinations(iKriging,inputIndices);
    ParacombinationsMatrix(end+1:end+size(Paracombinations,1),:) = Paracombinations;
    ambundanceVector(end+1:end+length(ambundance),1) = ambundance;
end
if obj.nMinimumDataPoint<0
    obj.ChosenCombinationsForPlot = ParacombinationsMatrix(ambundance==max(ambundance),:);
else
    obj.ChosenCombinationsForPlot = ambundanceVector(ambundance>obj.nMinimumDataPoint,:);
end

% Check Combinations
if isempty(obj.ChosenCombinationsForPlot)
    error('No appropriate combination of input variables are found. Please check class member "nMinimumDataPoint"')
end

%% Do Interpolation
doInterpolation(dimensionInterpolation);



%% Nested Functions
% -------------------------------------------------------------------------
function [] = checkLengthOfIndices(dimensionInterpolation)
    nIndicesDoesNotMatch = false;
    switch dimensionInterpolation
        case 2
            if length(inputIndices)~=1
                nIndicesDoesNotMatch = true;
            end
        case 3
            if length(inputIndices)~=2
                nIndicesDoesNotMatch = true;
            end
        otherwise
            error('dimensionInterpolation is only allowed to be 2 o 3 for 2D/3D interpolation,respectively')
    end
    if nIndicesDoesNotMatch
        error('length(indices) mus be %i',dimensionInterpolation-1)
    end
end
% -------------------------------------------------------------------------
function [] = doInterpolation(dimensionInterpolation)
    if isempty(indices2)
        % Make Sure no combination is chosen
        obj.ChosenCombinationsForPlot = [];
        switch dimensionInterpolation
            case 2
                obj.calcInterpolation_2D(KrigingObjectIndex,inputIndices);
                obj.plotInterpolation_2D(KrigingObjectIndex);   
            case 3
                obj.calcInterpolation_3D(KrigingObjectIndex,inputIndices);
                obj.plotInterpolation_3D(KrigingObjectIndex);   
            otherwise
                error('dimensionInterpolation is only allowed to be 2 o 3 for 2D/3D interpolation,respectively')
        end

    else
        switch dimensionInterpolation
            case 2
                for iComb = 1 : size(obj.ChosenCombinationsForPlot,1)
                    obj.calcInterpolation_2D(KrigingObjectIndex,inputIndices,indices2,obj.ChosenCombinationsForPlot(iComb,:));
                    obj.plotInterpolation_2D(KrigingObjectIndex);
                end
            case 3
                for iComb = 1 : size(obj.ChosenCombinationsForPlot,1)
                    obj.calcInterpolation_3D(KrigingObjectIndex,inputIndices,indices2,obj.ChosenCombinationsForPlot(iComb,:));
                    obj.plotInterpolation_3D(KrigingObjectIndex);
                end
            otherwise
                error('dimensionInterpolation is only allowed to be 2 o 3 for 2D/3D interpolation,respectively')
        end
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
