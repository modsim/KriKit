classdef TaylorKriging<KrigingSuperClass
% In this variant of Kriging the nonlinear parameter in the basis function 
% are estimated using a linear approaximation of the basis function.
% Only one basis function can be provided
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.
    
    
    %% Public Members
    properties(GetAccess='public',SetAccess='public')
    end
    
    %% Private Members
    properties(GetAccess='private',SetAccess='private')
        % Number of iterations to estimate the parameters in NonLinearKrigingLocalParaEstimation
        nIterationParaEstimation = 1;
%         nIterationParaEstimationInnerLoop = 1;
        % 
        BasisFctDerivative={};
        
        % The the order of the taylor expansion 
        TaylorExpansionOrder = 1;
        
        % Contains the intermidiate result during the parameter estimation
        % row represents iteration step -1
        SaveIntermediateResultsParaEstimation = [];
        
        % Decide if intermediate results should be saved
        SaveItermediateResults = 0;
        
        nBasisFctDerivative = 0;
        
        % Decide if the optimizer should use a auxiliary condition of
        % sum(Cinv(i,:).^2)>=ThresholdCondition
        UseAuxiliaryCondition = 0;
        % Determine the threshold in the auxiliary condition of
        % sum(Cinv(i,:).^2)>=ThresholdCondition
        ThresholdCondition = 0;
        
        
        InvCoVar = [];
        
        % Contains the Taylor coeffcientens estimated by Kriging using user
        % defined first derivative and a parameter set as linearization
        % point
        TaylorCoefficients = [];
        % Contains the Taylor coeffcientens estimated by Kriging using user
        % defined first derivative and a parameter set as linearization
        % point
        BasisFunctionCovarianceMatrix = [];
    end
    
    %% Protected Members
    properties(GetAccess='protected',SetAccess='protected')
    end
    
    %% Methods
    methods
        %% Initialization
        % Constructor
        function obj = TaylorKriging(varargin)
%             obj.EstimationVariogramType = 1;
%             obj.PredictionType = 1;
%             obj.MakePreworkType = 1;
        end
        

        %% ParameterEstimation
%         [] = estimateBasisFunctionParametersViaKriging(obj)
        [] = estimateBasisFctParametersViaTaylorKriging(obj)
        % -----------------------------------------------------------------
        [] = estimateBasisFctParametersFminCon(obj)
        % -----------------------------------------------------------------
        [] = estimateBasisFctParametersFminUnc(obj)
        % -----------------------------------------------------------------
        [] = estimateBasisFctParametersFminSearch(obj)
        % -----------------------------------------------------------------
        [] = estimateBasisFctParametersLocal(obj)
        % -----------------------------------------------------------------
        [] = estimateBasisFctParametersGA(obj)
        % -----------------------------------------------------------------
        [quality] = KrigingParaEstimationObjFct(obj,parameters)
        % -----------------------------------------------------------------
        [quality] = calcKrigingParaEstimationObjFctVector(obj,parameters)
        % -----------------------------------------------------------------
        function []=reset(obj)
            reset@KrigingSuperClass(obj);
%             obj.EstimationVariogramType = 1;
            obj.PredictionType = 4;
%             obj.MakePreworkType = 1;
            obj.nIterationParaEstimation = 1;
            obj.BasisFctDerivative={};
            obj.TaylorExpansionOrder = 1;
        end
        
        %% Get Function
        function derivative = getBasisFctDerivative(obj)
            derivative = obj.BasisFctDerivative;
        end
        %------------------------------------------------------------------
%         function TaylorExpansionOrder = getTaylorExpansionOrder(obj)
%             TaylorExpansionOrder = obj.TaylorExpansionOrder;
%         end
        
        % -----------------------------------------------------------------
        function [nIterationParaEstimation]=getnIterationParaEstimation(obj)
             nIterationParaEstimation = obj.nIterationParaEstimation ;
        end
        % -----------------------------------------------------------------
        function [SaveIntermediateResultsParaEstimation]=getSaveIntermediateResultsParaEstimation(obj)
            SaveIntermediateResultsParaEstimation = obj.SaveIntermediateResultsParaEstimation;
        end
        % -----------------------------------------------------------------
        function [SaveItermediateResults]=getSaveItermediateResults(obj)
            SaveItermediateResults = obj.SaveItermediateResults;
        end
        % -----------------------------------------------------------------
        function [TaylorExpansionOrder]=getTaylorExpansionOrder(obj)
            TaylorExpansionOrder = obj.TaylorExpansionOrder;
        end
        % -----------------------------------------------------------------
        function [nBasisFctDerivative]=getnBasisFctDerivative(obj)
            nBasisFctDerivative = obj.nBasisFctDerivative;
        end
        % -----------------------------------------------------------------
        function [UseAuxiliaryCondition]=getUseAuxiliaryCondition(obj)
            UseAuxiliaryCondition = obj.UseAuxiliaryCondition;
        end
        % -----------------------------------------------------------------
        function [ThresholdCondition]=getThresholdCondition(obj)
            ThresholdCondition = obj.ThresholdCondition;
        end
        % -----------------------------------------------------------------
        function [TaylorCoefficients]=getTaylorCoefficients(obj)
            TaylorCoefficients = obj.TaylorCoefficients;
        end
        % -----------------------------------------------------------------
        function [BasisFunctionCovarianceMatrix]=getBasisFunctionCovarianceMatrix(obj)
            BasisFunctionCovarianceMatrix = obj.BasisFunctionCovarianceMatrix;
        end
        % -----------------------------------------------------------------
        function [InvCoVar]=getInvCoVar(obj)
            InvCoVar = obj.InvCoVar;
        end
        
        %% Set Function
        function []=setBasisFctDerivative(obj,derivative)
            % setBasisFctDerivative(derivative)
            % "derivative" is a cell array containing the derivation of the
            % basis function (B) with respect to its parameters. Derivative
            % should be structured as following:
            % derivative={dB/dp_1;...;dB/dp_m;
            %              d^2B/(dp_1)^2;...;d^2B/(dp_m)^2;
            %              d^2B/(dp_1*dp_2);d^2B/(dp_1*dp_3);...;d^2B/(dp_1*dp_m);
            %              d^2B/(dp_2*dp_3);d^2B/(dp_1*dp_4);...;d^2B/(dp_2*dp_m);
            %              ...
            %              d^2B/(dp_(m-1)*dp_m)}
            % Note that the second derivatives start with the are at the
            % beginning twice derived with respect to the same parameter. After all
            % possible pairs are follwing.
            % Each derivative has to be formulate such as in the folllwoing example:
            % dB/(dp_1)='@(p,x)(x(:,1).*x(:,2).*x(:,5))./(p(2)+p(3)*x(:,1)+p(4)*x(:,2)+x(:,1).*x(:,2))'
            % You can set:
            % TaylorExpansionOrder ... degree of taylor polynomial which
            % may be created in other functions and use the derivatives
            
%             if size(derivative,1)~=obj.nBasisFctParameters&&size(derivative,2)==obj.nBasisFctParameters&&size(derivative,1)~=1
            if ischar(derivative)
                derivative = {derivative};
            end
            
            % Expect one derivative for each parameter (first order)
            % Expect one derivative for each parameter (second order)
            % Expect one derivative for each pair of parameter (second order)
            nDerivativeExpected = obj.nBasisFctParameters+(obj.nBasisFctParameters*(obj.nBasisFctParameters-1)/2 + obj.nBasisFctParameters)*(obj.TaylorExpansionOrder==2);
            
            if size(derivative,1)==1&&size(derivative,2)>=obj.nBasisFctParameters&&size(derivative,2)~=1
                warning('Input for setBasisFctDerivative is in a row vector, a column vector is expected. The provided vector is transposed');
                derivative = derivative';
            end
            if size(derivative,1)>=obj.nBasisFctParameters&&size(derivative,2)>1
                error('Input for setBasisFctDerivative should be of the dimension %dx1',obj.nBasisFctParameters);
%             elseif size(derivative,1)~=nDerivativeExpected&&size(derivative,2)==1&&obj.TaylorExpansionOrder==2
            elseif size(derivative,1)~=nDerivativeExpected&&size(derivative,2)==1
                warning('Expect excact %i but %i are provided (TaylorExpansionOrder = %i)',...
                    nDerivativeExpected,size(derivative,1),obj.TaylorExpansionOrder);
            elseif size(derivative,1)==nDerivativeExpected&&size(derivative,2)==1||obj.TaylorExpansionOrder==2&&size(derivative,1)~=obj.nBasisFctParameters*(obj.nBasisFctParameters-1)/2+obj.nBasisFctParameters
%                 warning('Fine')
                % Everything is fine
            else
                error('Something wrong in setBasisFctDerivative')
            end
            obj.BasisFctDerivative = cell(length(derivative),1);
            obj.nBasisFctDerivative = size(obj.BasisFctDerivative,1);
            
            for iPara = 1:length(derivative)
                try
                    obj.BasisFctDerivative{iPara}=eval(derivative{iPara});
                catch
                    warning('something wrong with your provided string: %s',derivative{iPara})
                    obj.BasisFctDerivative{iPara}=eval(derivative{iPara});
                end
            end
        end
        % -----------------------------------------------------------------
        function []=setnIterationParaEstimation(obj,nIterationParaEstimation)
            obj.nIterationParaEstimation = nIterationParaEstimation;
        end
        % -----------------------------------------------------------------
        function []=setSaveItermediateResults(obj,SaveItermediateResults)
            obj.SaveItermediateResults = SaveItermediateResults;
        end
        % -----------------------------------------------------------------
        function []=setTaylorExpansionOrder(obj,TaylorExpansionOrder)
            obj.TaylorExpansionOrder = TaylorExpansionOrder;
        end
        % -----------------------------------------------------------------
        function []=setUseAuxiliaryCondition(obj,UseAuxiliaryCondition)
            obj.UseAuxiliaryCondition = UseAuxiliaryCondition;
        end
        % -----------------------------------------------------------------
        function []=setThresholdCondition(obj,ThresholdCondition)
            obj.ThresholdCondition = ThresholdCondition;
        end
        % -----------------------------------------------------------------
        function []=setTaylorCoefficients(obj,TaylorCoefficients)
            obj.TaylorCoefficients = TaylorCoefficients;
        end
        % -----------------------------------------------------------------
        function []=setBasisFunctionCovarianceMatrix(obj,BasisFunctionCovarianceMatrix)
            obj.BasisFunctionCovarianceMatrix = BasisFunctionCovarianceMatrix;
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
