classdef BayesianOptimization<AnalyzeKriging
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.
    
    %% Public Members
    properties(GetAccess='public',SetAccess='public')
    end
    
    %% Private Members
    properties(GetAccess='private',SetAccess='private')
        ObjectiveIndicesUsedByCalcNewSamplesViaMCMC = [];
    end
    
    %% Protected Members
    properties(GetAccess='protected',SetAccess='protected')
        InequalityConstraintHandle = [];
        InequalityConstraintOutputHandle = [];
        nCutLinks = 0;
        ConsiderOnlyMaxExpectedImprovement = false;
    end
    
    methods
        %% Constructor
        function obj = BayesianOptimization()
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
        [newSamplePoints]=calcNewSamplesViaMCMC(obj,varargin)
        % ----------------------------------------------------------------
        [probabilityDensity]=MCMCDistributionFctDRAM(obj,varargin)
        % ----------------------------------------------------------------
        [probabilityDensity]=MCMCDistributionFctSlice(obj,varargin)
        %% Get Functions
        function [nMCMCLinks] = getnMCMCLinks(obj)
            nMCMCLinks = obj.nMCMCLinks;
        end
        % ----------------------------------------------------------------
        function [nCutLinks] = getnCutLinks(obj)
            nCutLinks = obj.nCutLinks;
        end
        % ----------------------------------------------------------------
        function [InequalityConstraintHandle] = getInequalityConstraintHandle(obj)
            InequalityConstraintHandle = obj.InequalityConstraintHandle;
        end
        % ----------------------------------------------------------------
        function [InequalityConstraintOutputHandle] = getInequalityConstraintOutputHandle(obj)
            InequalityConstraintOutputHandle = obj.InequalityConstraintOutputHandle;
        end
        % ----------------------------------------------------------------
        function [ConsiderOnlyMaxExpectedImprovement] = getConsiderOnlyMaxExpectedImprovement(obj)
            ConsiderOnlyMaxExpectedImprovement = obj.ConsiderOnlyMaxExpectedImprovement;
        end
        %% Set Functions
        function [] = setnMCMCLinks(obj,nMCMCLinks)
            if nMCMCLinks<=0||mod(nMCMCLinks,1)~=0
                error('nMCMCLinks has an positive integer bigger than 0')
            end
            obj.nMCMCLinks = nMCMCLinks;
        end
        % ----------------------------------------------------------------
        function [] = setInequalityConstraintHandle(obj,InequalityConstraintHandle)
            if isa(InequalityConstraintHandle,'function_handle')
                obj.InequalityConstraintHandle = InequalityConstraintHandle;
            else
                error('InequalityConstraintHandle has to be a function handle')
            end
        end
        % ----------------------------------------------------------------
        function [] = setInequalityConstraintOutputHandle(obj,InequalityConstraintOutputHandle)
            if isa(InequalityConstraintOutputHandle,'function_handle')
                obj.InequalityConstraintOutputHandle = InequalityConstraintOutputHandle;
            else
                error('InequalityConstraintHandle has to be a function handle')
            end
        end
        
        % ----------------------------------------------------------------
        function [] = setnCutLinks(obj,nCutLinks)
            if nCutLinks<=0||mod(nCutLinks,1)~=0
                error('nCutLinks has an positive integer bigger than 0')
            end
            obj.nCutLinks = nCutLinks;
        end
        % ----------------------------------------------------------------
        function [] = setConsiderOnlyMaxExpectedImprovement(obj,ConsiderOnlyMaxExpectedImprovement)
            if ~islogical(ConsiderOnlyMaxExpectedImprovement)
                error('ConsiderOnlyMaxExpectedImprovement has to be logical')
            end
            obj.ConsiderOnlyMaxExpectedImprovement = ConsiderOnlyMaxExpectedImprovement;
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
