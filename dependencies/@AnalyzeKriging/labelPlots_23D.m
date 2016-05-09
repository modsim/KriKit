function labelPlots_23D(obj,varargin)
% labelPlots_23D(obj,varargin)
% 
% label the 2D/3D interpolation plot. For Further details about the input
% see documentation of "plotInterpolation_23D()"
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

%% Initialization
KrigingObjectIndex = varargin{1};
dimensionInterpolation = varargin{2};
plotExpectedImprovement = varargin{3};
   
%% Actual Labeling
switch dimensionInterpolation
    case 2
        if isempty(obj.InputVarNames{KrigingObjectIndex(1)})
            xlabel(horzcat('Input Variable ',num2str(obj.KrigingPrediction_Interpolation2D{KrigingObjectIndex,3}(1))),'FontSize',20);
        else
            xlabel(obj.InputVarNames{KrigingObjectIndex(1)}(obj.KrigingPrediction_Interpolation2D{KrigingObjectIndex,3}(1)),'FontSize',20);
        end
        if plotExpectedImprovement
            ylabel(sprintf('%s\n%s','Exp. Improvement',obj.KrigingObjectNames{KrigingObjectIndex}),'FontSize',20);
        else
            ylabel(obj.KrigingObjectNames{KrigingObjectIndex},'FontSize',20);
        end
    case 3
%         title(obj.KrigingObjectNames{KrigingObjectIndex},'FontSize',obj.FontSize)

        if isempty(obj.InputVarNames{KrigingObjectIndex(1)})
            xlabel(horzcat('Input Variable ',num2str(obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,3}(1))),'FontSize',obj.FontSize);
            ylabel(horzcat('Input Variable ',num2str(obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,3}(2))),'FontSize',obj.FontSize);
        else
            xlabel(obj.InputVarNames{KrigingObjectIndex(1)}(obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,3}(1)),'FontSize',obj.FontSize);
            ylabel(obj.InputVarNames{KrigingObjectIndex(1)}(obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,3}(2)),'FontSize',obj.FontSize);
        end
        if plotExpectedImprovement
            zlabel('Expected Improvement','FontSize',obj.FontSize);
        else
            zlabel(obj.KrigingObjectNames{KrigingObjectIndex},'FontSize',obj.FontSize);
        end
    otherwise
        error('Plotting function is only allowed for 2D and 3D interpolation')
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
