function [inputData,indexValid]=plotOptimum23D(obj,varargin)
% [inputData,indexValid] = plotOptimum23D(obj,KrigingIndex,dimensionInterpolation,testValue)
% 
% Before this function can be used a 2D/3D interpolation has to be calculated!
% Use e.g. "calcInterpolation_2D"/"calcInterpolation_3D"
% Red stars indicate potentially optimum and blue stars indicate points
% which are significant suboptimal with a significance level of
% alpha="SignificanceLevel"
% This test ist based on a z-test. "inputData" contains sample points at which
% z-test was applied. "indexValid" is a bool vector, each entrie is
% associated with the associated entrie in inputData and is true if this
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
%     testValue = varargin{3};
    
    % In case of mutual Kriging interpolation, the prediction of all
    % objects is saved redundant, such that only the first index should
    % be used
    KrigingObjectIndex = obj.checkKrigingIndizes(KrigingObjectIndex);
    KrigingObjectIndexArray = KrigingObjectIndex;
    KrigingObjectIndex = KrigingObjectIndexArray(1);
    

    [inputData,indexValid]=doTestForOptimality(obj,varargin{:});
    indexInvalid = ~indexValid;
    
    names = obj.getInputVarNames(KrigingObjectIndex(1));
    switch dimensionInterpolation
        case 2
            prediction = obj.KrigingPrediction_Interpolation2D{KrigingObjectIndex,1}(:,1);
            InputVarIndices = obj.KrigingPrediction_Interpolation2D{KrigingObjectIndex,3};
            % Start with standard 2D Interpolation
            obj.plotInterpolation_2D(KrigingObjectIndex)
            hold on
            % Fill in the t-test results
            h_yes = plot(inputData(indexValid),prediction(indexValid),'r*');
            h_no = plot(inputData(indexInvalid),prediction(indexInvalid),'b*');
            % Name the axis
            if isempty(names)
                names = obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex(1),3};
                xlabel(num2str(names(InputVarIndices(1))),'FontSize',20)
            else
                xlabel(names{InputVarIndices(1)},'FontSize',20)
            end
            ylabel(obj.getKrigingObjectNames{KrigingObjectIndex},'FontSize',20);
        case 3
            InputVarIndices = obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,3};
            figure()
            hold on
            % Fill in the t-test results
            h_yes = plot(inputData(indexValid,1),inputData(indexValid,2),'r*');
            h_no = plot(inputData(indexInvalid,1),inputData(indexInvalid,2),'b*');
            % Name the axis
            if isempty(names)
                names = obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex(1),3};
                xlabel(num2str(names(InputVarIndices(1))),'FontSize',20)
                ylabel(num2str(names(InputVarIndices(2))),'FontSize',20)
            else
                xlabel(names{InputVarIndices(1)},'FontSize',20)
                ylabel(names{InputVarIndices(2)},'FontSize',20)
            end
            title(obj.getKrigingObjectNames{KrigingObjectIndex},'FontSize',20);
        otherwise
            error('Diminsion is unusual')
    end
    
    
    
    legend([h_yes,h_no],'Potential Optimum','No Optimum')

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
