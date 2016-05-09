classdef AnalyzeMixture<AnalyzeKriging
%
% For triangular plotting the Ternplot package was used available at 
% http://www.mathworks.com/matlabcentral/fileexchange/2299-ternploton
% Acces Date (02-02-2016)
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

	%% Public Members

    properties(GetAccess='public',SetAccess='public')
    end
    
    %% Private Members
    properties(GetAccess='private',SetAccess='private')
        
    end
    
    %% Protected Members
    properties(GetAccess='protected',SetAccess='protected')
    end
    
    methods
        %% Constructor
        function obj = AnalyzeMixture()
            obj@AnalyzeKriging();
        end
        %% Copy Operator for a shallow copy
        % ----------------------------------------------------------------
        function copy = copyObj(obj)
        % Create a shallow copy of the calling object.
            copy = eval(class(obj));
            meta = eval(['?',class(obj)]);
            for p = 1: size(meta.Properties,1)
                    pname = meta.Properties{p}.Name;
                try
                    eval(['copy.',pname,' = obj.',pname,';']);
                catch
                    error(['\nCould not copy ',pname,'.\n']);
%                     fprintf(['\nCould not copy ',pname,'.\n']);
                end
            end
        end
        %% General Methods
        % ----------------------------------------------------------------
        function [inputData,indexValid]=plotOptimum3D(obj,varargin)
            fctName = dbstack('-completenames');
            error('%s does not exist',fctName.name)
        end
        % ----------------------------------------------------------------
        function [inputData,indexValid]=plotOptimum2D(obj,varargin)
            fctName = dbstack('-completenames');
            error('%s does not exist',fctName.name)
        end
        % ----------------------------------------------------------------
        function  [inputData,indexValid]=plotOptimum23D(obj,varargin)
            fctName = dbstack('-completenames');
            error('%s does not exist',fctName.name)
        end
        % ----------------------------------------------------------------
        function []=calcAndPlotInterpolation_2D_BestChoice(obj,varargin)
            fctName = dbstack('-completenames');
            error('%s does not exist',fctName.name)
        end
        % ----------------------------------------------------------------
        function []=calcInterpolation_3D(obj,varargin)
            fctName = dbstack('-completenames');
            error('%s does not exist',fctName.name)
        end
        % ----------------------------------------------------------------
        function [expectedParetoCurve,deviationParetoCurve,globalUncertainity] = predictParetoCurve(obj,varargin)
            fctName = dbstack('-completenames');
            error('%s does not exist',fctName.name)
        end
        % ----------------------------------------------------------------
        function []=calcParetoFront(obj,varargin)
            fctName = dbstack('-completenames');
            error('%s does not exist',fctName.name)
        end
        % ----------------------------------------------------------------
        function []=plotParetoFront(obj,varargin)
            fctName = dbstack('-completenames');
            error('%s does not exist',fctName.name)
        end
        % ----------------------------------------------------------------
        function []=plotRatio(obj,varargin)
            fctName = dbstack('-completenames');
            error('%s does not exist',fctName.name)
        end
        % ----------------------------------------------------------------
        function []=plotParetoInput_2D(obj,varargin)
            fctName = dbstack('-completenames');
            error('%s does not exist',fctName.name)
        end
        % ----------------------------------------------------------------
        function [] = calcnNewSamplesParetoOptimal(obj,varargin)
            fctName = dbstack('-completenames');
            error('%s does not exist',fctName.name)
        end
        % ----------------------------------------------------------------
        function [] = calcnNewSamples(obj,varargin)
            fctName = dbstack('-completenames');
            error('%s does not exist',fctName.name)
        end
        % ----------------------------------------------------------------
        function []=calcInterpolation_2D(obj,varargin)
            error('calcInterpolation_2D does not exisit in class AnalyzeMixture')
        end
        % ----------------------------------------------------------------
        function [] = calcnNewSamplesMarkovChain(obj,varargin)
            fctName = dbstack('-completenames');
            error('%s does not exist',fctName.name)
        end
        % ----------------------------------------------------------------
        function [ExpectedImprovement] = calcExpectedImprovementPareto(obj,varargin)
            fctName = dbstack('-completenames');
            error('%s does not exist',fctName.name)
        end
        % ----------------------------------------------------------------
        function [] = make3dMovieAnalysis(obj,KrigingObjectIndex)
            fctName = dbstack('-completenames');
            error('%s does not exist',fctName.name)
        end
        % ----------------------------------------------------------------
        function [] = calcScreeningAnalysis(obj,KrigingObjectIndex)
            fctName = dbstack('-completenames');
            error('%s does not exist',fctName.name)
        end
        % ----------------------------------------------------------------
        function [] = plotScreeningAnalysis(obj,varargin)
            fctName = dbstack('-completenames');
            error('%s does not exist',fctName.name)
        end
        % ----------------------------------------------------------------
        function [] = plotScreeningAnalysisKrigingInterpolation(obj,varargin)
            fctName = dbstack('-completenames');
            error('%s does not exist',fctName.name)
        end
        % ----------------------------------------------------------------
        function [] = plotScreeningAnalysisExpectedImprovement(obj,varargin)
            fctName = dbstack('-completenames');
            error('%s does not exist',fctName.name)
        end
        % ----------------------------------------------------------------
        []=calcInterpolation(obj,varargin);
        % ----------------------------------------------------------------
        [] = plotInterpolation_nD(obj,varargin);
        % ----------------------------------------------------------------
        function []=plotInterpolation_3D(obj,KrigingObjectIndex)
            % []=plotInterpolation_3D(KrigingObjectIndex)
            %
            % For further Details see documentation of "plotInterpolation_23D()"
            %
            % For triangular plotting the Ternplot package was used available at 
            % http://www.mathworks.com/matlabcentral/fileexchange/2299-ternploton
            % Acces Date (02-02-2016)
            obj.plotInterpolation_23D(KrigingObjectIndex,3);
        end
        % ----------------------------------------------------------------
        function []=plotInterpolation_2D(obj,KrigingObjectIndex)
            % []=plotInterpolation_2D(KrigingObjectIndex)
            %
            % For further Details see documentation of "plotInterpolation_23D()"
            %
            % For triangular plotting the Ternplot package was used available at 
            % http://www.mathworks.com/matlabcentral/fileexchange/2299-ternploton
            % Acces Date (02-02-2016)
            obj.plotInterpolation_23D(KrigingObjectIndex,2);
        end
        % ----------------------------------------------------------------
        [Data,dataShown,OutlierShown]=determineRelevantDataPointsForPlot(obj,varargin)
        % ----------------------------------------------------------------
        [KrigingObjectIndex,InputVar1,InputVar2,InputVar3,RemainingIndices,RemainingValues] = checkInputInterpolation3D(obj,varargin);
        % ----------------------------------------------------------------
        [] = labelPlots_23D(obj,varargin);
        % ----------------------------------------------------------------
        [] = plotInterpolation_23D(obj,varargin)
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
