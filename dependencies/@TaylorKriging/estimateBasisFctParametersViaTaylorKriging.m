function [] = estimateBasisFctParametersViaTaylorKriging(obj)
% Estimation of basis function is done iteratively. Each round the a
% taylor expansion is done. Following the basis function is
% approximated as sum of its derivatives and the coefficients contain
% the unknown parameters. Consequently the parameters are calculated
% out of the estimated taylor coefficients. However, snce taylor
% expansion is approximation only around the given starting parameter
% values this approach is redone several times using the parameter
% values which are estimated in the previous round
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

% ###### Initialization ###### 
if isempty(obj.InitialBasisFctParameters)
    warning('No initial values for the basis parameters are provided but parameter estimation via Kriging need initial values. Every initial value is set to one! You can change this using "setInitialBasisFctParameters"')
    obj.InitialBasisFctParameters = ones(obj.nBasisFctParameters,1);
end
obj.BasisFctParameters = obj.InitialBasisFctParameters;

if obj.ShowDetails==1
    InitialDisplay
end
% Save Results
if obj.SaveItermediateResults==1
    obj.SaveIntermediateResultsParaEstimation = zeros(obj.getnIterationParaEstimation+1,obj.getnBasisFctParameters+1);
    obj.estimateBasisFctCoefficients();
    obj.SaveIntermediateResultsParaEstimation(1,1)=obj.getBasisFctCoefficients;
    obj.SaveIntermediateResultsParaEstimation(1,2:end)=obj.BasisFctParameters';
end

% ###### Do checks ###### 
    % Check if enough expansion terms exists
if length(obj.getBasisFctDerivative)<obj.nBasisFctParameters
    error('Not enough derivation terms are provided. Derivative until the order %i are expected.',obj.getTaylorExpansionOrder*obj.nBasisFctParameters)
end
    % Check if more than one basis fct is defined
if obj.nBasisFct>1
    error('"estimateBasisFctParametersViaKriging" is only defined for one basis function which has only nonlinear parameters but %i basis function are defined',obj.nBasisFct)
end
    % Check if basis fct and its derivatives are linear independent from
    % each other
    indicesOfLinearDependentFcts = getIndicesOfLinearDependentFcts();
    if ~isempty(indicesOfLinearDependentFcts)
        indicesOfLinearDependentFcts-1
        warning('The above mentioned indices represent the derivatives which have the risk to be linear dependent to each other since their function values are very similar to each other for all input data values (0 means the basis function itself). Taylor Kriging does not work anymore if two or more functions are linear dependent to each other')
    end
    
    

    
% ###### Iterative estimation ###### 
for iIter= 1 : obj.getnIterationParaEstimation
    
    % Save Results
    if obj.SaveItermediateResults==1
        obj.estimateBasisFctCoefficients();
        obj.SaveIntermediateResultsParaEstimation(iIter+1,1)=obj.getBasisFctCoefficients;
    end
    
        % Calculate the coeffecicients of the taylor expansion
    obj.estimateBasisFctCoefficients(1);
    
    
        % Recalculate the parameters out of the coefficients
    for iPara = 1 : obj.nBasisFctParameters
        obj.BasisFctParameters(iPara)=obj.BasisFctCoefficients(iPara+1)/obj.BasisFctCoefficients(1)+obj.BasisFctParameters(iPara);
    end

        % Display Results
    if obj.ShowDetails==1
        IntermediateDisplay
    end
    
    % Save Results
    if obj.SaveItermediateResults==1
        obj.SaveIntermediateResultsParaEstimation(iIter+1,2:end)=obj.BasisFctParameters';
    end
end

% During "estimateBasisFctCoefficients" the extended covairance
% maxtrix was calculated this should be redone for the final parameters
obj.checkVariogram = 0;
%% ------------------------------------------------------------------------
    function [indicesOfLinearDependentFcts] = getIndicesOfLinearDependentFcts()
        % Initialization
            % Initialy, all fcts have the potential to be linear dependent
            % to each other
        nFcts = 1+obj.nBasisFctParameters;
        riskyIndices = 1:nFcts;
        
        % Collect the values of all fcts at the input data
        fctValues = zeros(obj.nExperiments,nFcts);
        for iFct = 1 : nFcts
            if iFct==1
                fctValues(:,1) = obj.BasisFct{1}(obj.BasisFctParameters,obj.InputData);
            else
                fctValues(:,iFct) = obj.BasisFctDerivative{iFct-1}(obj.BasisFctParameters,obj.InputData);
            end
        end
        
        % Test if no function is just a constant and show error if more
        % then one is constant ove all values, since then they are linear
        % dependent
            % Calculate the Covariance
        covMatrix = cov(fctValues);
            % Find zero entries = constant in all entries
        rConstant = find(abs(diag(covMatrix))<=1e-2);
            % Throw error if more then one constant fct
        if length(rConstant)>1
            rConstant-1
            indicesOfLinearDependentFcts=[];
            warning('The above mentioned indices represent the derivatives which are (almost) constant at all input data values (0 means the basis function itself). If more than 1 function show this behavior then theses function are linear dependent and the taylor Kriging does not work anymore')
            return
        else
            % Constant fct has no risk to be linear dependent with
            % function which values differ for differen input values
           riskyIndices= setdiff(riskyIndices,rConstant);
        end
        fctValues = fctValues(:,riskyIndices);
        
        % Find linear dependent functions by analysis its correlation to
        % each other
        corrMatrix = corrcoef(fctValues);
            % Set all irrelevant entrie equal to zero
        corrMatrix = triu(corrMatrix)- eye(length(riskyIndices));
            % Find all functions which have the risk to be linear dependent
            % to each other
        [r,c] = find(corrMatrix>0.98);
            % Since one column was erased the indices may shifted and have
            % to be adjusted
        if ~isempty(rConstant)
            r(r>=rConstant)=r(r>=rConstant)+1;
            c(c>=rConstant)=c(c>=rConstant)+1;
        end
        % Output 
        indicesOfLinearDependentFcts=[r,c];
        
    end
%% ------------------------------------------------------------------------
% Initial Display
    function []=InitialDisplay()
        fprintf('Iteration\t|Parameters\t|Estimation\n')
        for iParam = 1 : obj.nBasisFctParameters
            fprintf('%i\t|%i\t|%d\n',0,iParam,obj.BasisFctParameters(iParam));
        end
    end
%% ------------------------------------------------------------------------
function []=IntermediateDisplay()
    fprintf('-----------------------------------\n')
    fprintf('Constant Multiplicator: %d\n',obj.BasisFctCoefficients(1));
    for iParam = 1 : obj.nBasisFctParameters
        fprintf('%i\t|%i\t|%d\n',iIter,iParam,obj.BasisFctParameters(iParam));
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
