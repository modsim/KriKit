classdef OrdinaryKriging<KrigingSuperClass
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.
    %KRIGINGTOOL Summary of this class goes here
    %   Detailed explanation goes here
    
    %% Public Members
    properties(GetAccess='public',SetAccess='public')
    end
    
    %% Private Members
    properties(GetAccess='private',SetAccess='private')
    end
    
    %% Protected Members
    properties(GetAccess='protected',SetAccess='protected')
    end
    
    %% Methods
    methods
        %% Initialization
        % Constructor
        function obj = OrdinaryKriging(varargin)
            obj.BasisFct = cell(1,1);
            obj.BasisFct{1} = @(p,x)ones(size(x,1),1);
            obj.nBasisFct = 1;
            obj.BasisFctType = 'polynomial';
            obj.MaxDegree = 0;
%             obj.EstimationVariogramType = 1;
%             obj.PredictionType = 1;
%             obj.MakePreworkType = 1;
        end
        %% Calculate the inverse of the (Co)-Variogram
        %--------------------------------------
        % Set everything to the default values
%         reset(obj);
        %--------------------------------------
        calcInverseCoVar(obj);
        %% Model Prediction 
        % -----------------------------------------------------------------
        [output]=predictionMaxLikelihood(obj,input);
        %% Get Functions
        %% Set Functions
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
