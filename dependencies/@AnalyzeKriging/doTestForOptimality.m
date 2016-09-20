function [inputData,indexValid]=doTestForOptimality(obj,varargin)
% [inputData,indexValid] = doTestForOptimality(obj,KrigingIndex,dimensionInterpolation,testValue)
% 
% Before this function can be used a 2D/3D interpolation has to be calculated!
% Use e.g. "calcInterpolation_2D"/"calcInterpolation_3D"
% Red stars indicate potentially optimum and blue stars indicate points
% which are significant suboptimal with a significance level of
% alpha="SignificanceLevel"
% This test is based on a z-test. "inputData" contains sample points at
% which z-test was applied. "indexValid" is a bool vector, each entry is
% associated with the associated entry in inputData and is true if this
% point is accepted as part of the optimal region.
% 
% Input:
% KrigingIndex ... index of Kriging object of interest
% dimensionInterpolation ... decided if 2D(dimensionInterpolation=2) or
%                            3D(dimensionInterpolation=3) interpolation is
%                            considered
% testValue ... if provided, the z-test is done w.r.t to testValue
%
% Output:
% inputData ... contains sample points at which z-test was applied.
%               (nValidationPointsXnInputVar) 
% indexValid ... bool vector, each entry is associated with its
%                correspodent in inputData. If true, point is accepted as
%                part of the optimal region.
%
% You can set 
% - MinMax ... Vector which contains for each kriging objective 1 or -1
%              when the optimization goal is maximization or minimization,
%              respectively. By default minimization
% - useDataPointsAsComparisonPoint ... Decide if z-test shall be done
%                                      w.r.t best sample value. If false,
%                                      z-test is done w.r.t best
%                                      interpolation value (default false)
% - SignificanceLevel ... Signifcance level (error Erro Type I) used for z-test
%
% You can get : -
% 
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

    KrigingObjectIndex = varargin{1};
    dimensionInterpolation = varargin{2};
    testValue = varargin{3};
    
    % Decide which pair wise combination is considered in case of Screening
    if length(varargin)>3
        iCombination = varargin{4};
    else
        iCombination = 0;
    end
    
    % In case of mutual Kriging interpolation, the prediction of all
    % objects is saved redundant, such that only the first index should
    % be used
    KrigingObjectIndex = obj.checkKrigingIndizes(KrigingObjectIndex);
    KrigingObjectIndexArray = KrigingObjectIndex;
    KrigingObjectIndex = KrigingObjectIndexArray(1);
    

    % Get data and rescale if neccesary
    switch dimensionInterpolation
        case 2
%             InputVarIndices = obj.KrigingPrediction_Interpolation2D{KrigingObjectIndex,3};
            prediction = obj.KrigingPrediction_Interpolation2D{KrigingObjectIndex,1}(:,1);
            sd = obj.KrigingPrediction_Interpolation2D{KrigingObjectIndex,1}(:,2);
            inputVar1=(obj.KrigingPrediction_Interpolation2D{KrigingObjectIndex,2}(:,obj.KrigingPrediction_Interpolation2D{KrigingObjectIndex,3}(1)));
            inputData = inputVar1;
        case 3
%             InputVarIndices = obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,3};
            prediction = obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,1}(:,1);
            sd = obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,1}(:,2);
            inputVar1=(obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,2}(:,obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,3}(1)));
            inputVar2=(obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,2}(:,obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,3}(2)));
            inputData = [inputVar1,inputVar2];
        case 4 
            % Case Screening 
            prediction = obj.KrigingPrediction_Screening{KrigingObjectIndex,1}{iCombination}(:,1);
            sd = obj.KrigingPrediction_Screening{KrigingObjectIndex,1}{iCombination}(:,2);
            inputVarIndices = obj.KrigingPrediction_Screening{KrigingObjectIndex,3}{iCombination};
            inputVar1=(obj.KrigingPrediction_Screening{KrigingObjectIndex,2}{iCombination}(:,inputVarIndices(1)));
            inputVar2=(obj.KrigingPrediction_Screening{KrigingObjectIndex,2}{iCombination}(:,inputVarIndices(2)));
            inputData = [inputVar1,inputVar2];
        case 5
            % Case Screening 
            prediction = obj.KrigingPrediction_InterpolationnD{KrigingObjectIndex,1}(:,1);
            sd = obj.KrigingPrediction_InterpolationnD{KrigingObjectIndex,1}(:,2);
            inputData = obj.KrigingPrediction_InterpolationnD{KrigingObjectIndex,4};
        otherwise
            error('Diminsion is unusual')
    end

    % Do t-test with high degree of freedom
    minOrMax = obj.getMinMax(KrigingObjectIndex);
    if obj.UseDataPointsAsComparisonPoint
        [~,dataShown,~]=obj.determineRelevantDataPointsForPlot(KrigingObjectIndex,dimensionInterpolation);
        OutputData = obj.KrigingObjects{KrigingObjectIndex}.getOutputData;
        if sum(dataShown)==0
            error('If UseDataPointsAsComparisonPoint==true: Use only values for remaining input variables where at least one data point exist')
        end
    end
    
    if isempty(testValue)
        if minOrMax<0
            if obj.UseDataPointsAsComparisonPoint
                bestValue = min(OutputData(dataShown));
            else
                bestValue = min(prediction(:,1));
            end
        else
            if obj.UseDataPointsAsComparisonPoint
                bestValue = max(OutputData(dataShown));
            else
                bestValue = max(prediction(:,1));
            end
        end
    else
        bestValue = testValue;
    end

    
    z_score = (prediction(:,1)-bestValue)./sd;
    if minOrMax<0
        thresholdQuantile = sqrt(2)*erfinv(2*(1-obj.SignificanceLevel)-1);
        appropriatePoints = z_score>thresholdQuantile;
    else
        thresholdQuantile = sqrt(2)*erfinv(2*(obj.SignificanceLevel)-1);
        appropriatePoints = z_score<thresholdQuantile;
    end
    pTest = ones(size(prediction,1),1);
    pTest(appropriatePoints) = 0;
    
    % Plot figure
    indexValid = pTest==1;

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
