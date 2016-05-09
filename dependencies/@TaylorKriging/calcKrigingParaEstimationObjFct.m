function [quality] = calcKrigingParaEstimationObjFct(obj,para)
% Calculates the quality of the parameter estimation with respect to the
% Kriging quality measurement: log10(sum(a_l^2))
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

    % Create extended Covariance Matrix
        % Initialization
    if isempty(obj.CovariogramMatrix)
        obj.calcCovariogramMatrix
    end
    quality = 0;

    if length(para)~=obj.getnBasisFctParameters
        warning('Number of given parameters (%i) and the needed number of parameters (%i) differes from each other',length(para),obj.getnBasisFctParameters)
    end
    if obj.nBasisFctDerivative<(obj.getnBasisFctParameters+(obj.TaylorExpansionOrder>=2)*(obj.getnBasisFctParameters+obj.getnBasisFctParameters*(obj.getnBasisFctParameters-1)/2))
        error('Not enough derivatives are provided')
    end
    C = obj.CovariogramMatrix(1:obj.getnExperiments,1:obj.getnExperiments);
    F = zeros(obj.getnExperiments,1+obj.nBasisFctDerivative);
    
%     C(1:obj.getnExperiments,1:obj.getnExperiments) = ...
%                         obj.CovariogramMatrix(1:obj.getnExperiments,1:obj.getnExperiments);

        % Basis Function
    F(:,1) = obj.getBasisFct{1}(para,obj.getInputData);
        % Derivatives
    for iDerivative=1:obj.nBasisFctDerivative
        F(:,iDerivative+1) = obj.getBasisFctDerivative{iDerivative}(para,obj.getInputData);
    end
    
%     C=triu(C)+triu(C)'-bsxfun(@times,eye(length(C)),diag(C));
%     Cinv = inv(C);
    if isempty(obj.getInvCoVar)
        obj.InvCoVar = inv(C);
    end
%     saveCoeff = (F'*(C\F))\(F'*(C\(obj.getOutputData)));
    obj.BasisFunctionCovarianceMatrix = inv((F'*obj.InvCoVar*F));
%     saveCoeff = (F'*(obj.InvCoVar*F))\(F'*(obj.InvCoVar*(obj.getOutputData)));
    saveCoeff = obj.BasisFunctionCovarianceMatrix*(F'*(obj.InvCoVar*(obj.getOutputData)));
    quality = sum(saveCoeff(2:end).^2);
    quality = log10(quality);
    if isinf(quality)
        quality=-1e6;
    end
    

      obj.TaylorCoefficients = saveCoeff(2:end);
      
%     quality = sum(saveCoeff(2:end));
%     p0 = saveCoeff(1);
    
    
%     %% Check
%     % ####### Single Parameter estimation #######
%     % Create extended Covariance Matrix
%             % Initialization
%     if isempty(obj.CovariogramMatrix)
%         obj.calcCovariogramMatrix
%     end
%     quality = 0;
% 
%     if length(para)~=obj.getnBasisFctParameters
%         warning('Number of given parameters (%i) and the needed number of parameters (%i) differes from each other',length(para),obj.getnBasisFctParameters)
%     end
%     if obj.nBasisFctDerivative<(obj.getnBasisFctParameters+(obj.TaylorExpansionOrder>=2)*(obj.getnBasisFctParameters+obj.getnBasisFctParameters*(obj.getnBasisFctParameters-1)/2))
%         error('Not enough derivatives are provided')
%     end
%     C = zeros(obj.getnExperiments+1+obj.nBasisFctDerivative);
%     C(1:obj.getnExperiments,1:obj.getnExperiments) = ...
%                         obj.CovariogramMatrix(1:obj.getnExperiments,1:obj.getnExperiments);
% 
%         % Basis Function
%     C(1:obj.getnExperiments,obj.getnExperiments+1) = obj.getBasisFct{1}(para,obj.getInputData);
%         % Derivatives
%     for iDerivative=1:obj.nBasisFctDerivative
%         C(1:obj.getnExperiments,obj.getnExperiments+iDerivative+1) = obj.getBasisFctDerivative{iDerivative}(para,obj.getInputData);
%     end
% 
%     C=triu(C)+triu(C)'-bsxfun(@times,eye(length(C)),diag(C));
%     Cinv = inv(C);
%     saveCoeff = Cinv(obj.getnExperiments+2:end,:)*[obj.getOutputData;zeros(obj.nBasisFctDerivative+1,1)];
%     quality = sum(saveCoeff.^2);
%     quality = log10(quality);
end


%    C = zeros(obj.getnExperiments+1+obj.nBasisFctDerivative);
%     C(1:obj.getnExperiments,1:obj.getnExperiments) = ...
%                         obj.CovariogramMatrix(1:obj.getnExperiments,1:obj.getnExperiments);
% 
%         % Basis Function
%     C(1:obj.getnExperiments,obj.getnExperiments+1) = obj.getBasisFct{1}(para,obj.getInputData);
%         % Derivatives
%     for iDerivative=1:obj.nBasisFctDerivative
%         C(1:obj.getnExperiments,obj.getnExperiments+iDerivative+1) = obj.getBasisFctDerivative{iDerivative}(para,obj.getInputData);
%     end
%     
%         % All together
%     C=triu(C)+triu(C)'-bsxfun(@times,eye(length(C)),diag(C));
%         % Calculate the objective function
%         Cinv = inv(C);
% %         saveCoeff(iDerivative) =(Cinv(end,:)*[obj.getOutputData;zeros(2,1)])/sum(Cinv(end,:).^2);
%     for iDerivative=1:obj.nBasisFctDerivative
%         saveCoeff(iDerivative) =(Cinv(obj.getnExperiments+1+iDerivative,:)*[obj.getOutputData;zeros(1+obj.nBasisFctDerivative,1)]);
%     end

    
% % Calculate the objective function
%     Cinv = inv(C);
% 
% %     for iPara = 1:length(para)
% %         quality = quality + log10(sum((Cinv(obj.getnExperiments+iPara+1:end,:)*[obj.getOutputData;zeros(obj.getnBasisFctParameters+1,1)]*100).^2));
% %     end
% %     saveCoeff = zeros(obj.nBasisFctDerivative,1);
% %     for iDerivative=1:obj.nBasisFctDerivative
% %         saveCoeff(iDerivative) = sum((Cinv(obj.getnExperiments+iDerivative+1:end,:)*[obj.getOutputData;zeros(obj.nBasisFctDerivative+1,1)]*100).^2);
% %         quality = quality + log10(sum((Cinv(obj.getnExperiments+iDerivative+1:end,:)*[obj.getOutputData;zeros(obj.nBasisFctDerivative+1,1)]*100).^2));
% %     end


% ######################################################################
% ####### Single Parameter estimation #######
% % Create extended Covariance Matrix
%         % Initialization
%     if isempty(obj.CovariogramMatrix)
%         obj.calcCovariogramMatrix
%     end
%     quality = 0;
% 
%     if length(para)~=obj.getnBasisFctParameters
%         warning('Number of given parameters (%i) and the needed number of parameters (%i) differes from each other',length(para),obj.getnBasisFctParameters)
%     end
%     if obj.nBasisFctDerivative<(obj.getnBasisFctParameters+(obj.TaylorExpansionOrder>=2)*(obj.getnBasisFctParameters+obj.getnBasisFctParameters*(obj.getnBasisFctParameters-1)/2))
%         error('Not enough derivatives are provided')
%     end
%     C = zeros(obj.getnExperiments+1+obj.nBasisFctDerivative);
%     C(1:obj.getnExperiments,1:obj.getnExperiments) = ...
%                         obj.CovariogramMatrix(1:obj.getnExperiments,1:obj.getnExperiments);
% 
%         % Basis Function
%     C(1:obj.getnExperiments,obj.getnExperiments+1) = obj.getBasisFct{1}(para,obj.getInputData);
%         % Derivatives
%     for iDerivative=1:obj.nBasisFctDerivative
%         C(1:obj.getnExperiments,obj.getnExperiments+iDerivative+1) = obj.getBasisFctDerivative{iDerivative}(para,obj.getInputData);
%     end
%     
%     C=triu(C)+triu(C)'-bsxfun(@times,eye(length(C)),diag(C));
%     Cinv = inv(C);
%     saveCoeff = Cinv(obj.getnExperiments+2:end,:)*[obj.getOutputData;zeros(obj.nBasisFctDerivative+1,1)];
%     quality = sum(saveCoeff.^2);
%     quality = log10(quality);



%     for iDerivative=1:obj.nBasisFctDerivative
%         C = zeros(obj.getnExperiments+1+1);
%         C(1:obj.getnExperiments,1:obj.getnExperiments) = ...
%                             obj.CovariogramMatrix(1:obj.getnExperiments,1:obj.getnExperiments);
% 
%             % Basis Function
%         C(1:obj.getnExperiments,obj.getnExperiments+1) = obj.getBasisFct{1}(para,obj.getInputData);
%             % Derivatives
%             C(1:obj.getnExperiments,end) = obj.getBasisFctDerivative{iDerivative}(para,obj.getInputData);
%             % All together
%         C=triu(C)+triu(C)'-bsxfun(@times,eye(length(C)),diag(C));
% 
%     % Calculate the objective function
%         Cinv = inv(C);
%         saveCoeff(iDerivative) =(Cinv(end,:)*[obj.getOutputData;zeros(2,1)]);
%     end
%       quality = log10(sum(saveCoeff.^2));


%     X = zeros(obj.getnExperiments,1+obj.nBasisFctDerivative);
%     X(:,1) = obj.getBasisFct{1}(para,obj.getInputData);
%     
%     C = zeros(obj.getnExperiments+1+obj.nBasisFctDerivative);
%     C(1:obj.getnExperiments,1:obj.getnExperiments) = ...
%                         obj.CovariogramMatrix(1:obj.getnExperiments,1:obj.getnExperiments);
% 
%         % Basis Function
%     C(1:obj.getnExperiments,obj.getnExperiments+1) = obj.getBasisFct{1}(para,obj.getInputData);
%         % Derivatives
%     for iDerivative=1:obj.nBasisFctDerivative
%         C(1:obj.getnExperiments,obj.getnExperiments+iDerivative+1) = obj.getBasisFctDerivative{iDerivative}(para,obj.getInputData);
%     end
%         % All together
%     C=triu(C)+triu(C)'-bsxfun(@times,eye(length(C)),diag(C));
% 
% % Calculate the objective function
%     Cinv = inv(C);
% 
% %     for iPara = 1:length(para)
% %         quality = quality + log10(sum((Cinv(obj.getnExperiments+iPara+1:end,:)*[obj.getOutputData;zeros(obj.getnBasisFctParameters+1,1)]*100).^2));
% %     end
% % %     saveCoeff = zeros(obj.nBasisFctDerivative,1);
%     for iDerivative=1:obj.nBasisFctDerivative
%         saveCoeff(iDerivative) = (Cinv(obj.getnExperiments+iDerivative+1,:)*[obj.getOutputData;zeros(obj.nBasisFctDerivative+1,1)]);
%     end

%%
% if obj.nBasisFctParameters>obj.TaylorExpansionOrder
%     comb = nchoosek(1:obj.nBasisFctParameters,obj.TaylorExpansionOrder);
% else
%     comb = 1:obj.nBasisFctParameters;
% end
% ncomb = size(comb,1);
% % Initialize covariance matrix: covariance matrix + basis function
% % plus first derivatives with respect to the two parameters and its
% % combination (df^2/dp1/dp2)
% 
% 
% if obj.TaylorExpansionOrder==2
%     saveCoeff = zeros(ncomb*3-2*obj.nBasisFctParameters,1);
%     C = zeros(obj.getnExperiments+1+3);
%     C(1:obj.getnExperiments,1:obj.getnExperiments) = ...
%                             obj.CovariogramMatrix(1:obj.getnExperiments,1:obj.getnExperiments);
%     for iComb=1:ncomb
% 
%         % Basis Function
%         C(1:obj.getnExperiments,obj.getnExperiments+1) = obj.getBasisFct{1}(para,obj.getInputData);
%         % Derivatives
%         C(1:obj.getnExperiments,obj.getnExperiments+2) = obj.getBasisFctDerivative{comb(iComb,1)}(para,obj.getInputData);
%         C(1:obj.getnExperiments,obj.getnExperiments+3) = obj.getBasisFctDerivative{comb(iComb,2)}(para,obj.getInputData);
%         C(1:obj.getnExperiments,obj.getnExperiments+4) = obj.getBasisFctDerivative{obj.nBasisFctParameters+iComb}(para,obj.getInputData);
%         % Make it symmectrical
%         C(obj.getnExperiments+1,1:obj.getnExperiments) = obj.getBasisFct{1}(para,obj.getInputData)';
%         C(obj.getnExperiments+2,1:obj.getnExperiments) = obj.getBasisFctDerivative{comb(iComb,1)}(para,obj.getInputData)';
%         C(obj.getnExperiments+3,1:obj.getnExperiments) = obj.getBasisFctDerivative{comb(iComb,2)}(para,obj.getInputData)';
%         C(obj.getnExperiments+4,1:obj.getnExperiments) = obj.getBasisFctDerivative{obj.nBasisFctParameters+iComb}(para,obj.getInputData)';
% 
%         % All together
% %         C=triu(C)+triu(C)'-bsxfun(@times,eye(length(C)),diag(C));
% 
%         % Calculate the objective function
%         Cinv = inv(C);
%         saveCoeff((iComb-1)*3+1) =(Cinv(obj.getnExperiments+2,:)*[obj.getOutputData;zeros(2,1)]);
%         saveCoeff((iComb-1)*3+2) =(Cinv(obj.getnExperiments+3,:)*[obj.getOutputData;zeros(2,1)]);
%         saveCoeff((iComb-1)*3+3) =(Cinv(obj.getnExperiments+4,:)*[obj.getOutputData;zeros(2,1)]);
%     end
% else
%     C = zeros(obj.getnExperiments+1+1);
%     C(1:obj.getnExperiments,1:obj.getnExperiments) = ...
%                             obj.CovariogramMatrix(1:obj.getnExperiments,1:obj.getnExperiments);
%     for iPara=1:obj.BasisFctParameters
%         % Basis Function
%         C(1:obj.getnExperiments,obj.getnExperiments+1) = obj.getBasisFct{1}(para,obj.getInputData);
%         % Derivatives
%         C(1:obj.getnExperiments,obj.getnExperiments+2) = obj.getBasisFctDerivative{comb(iPara,1)}(para,obj.getInputData);
%         % Make it symmectrical
%         C(obj.getnExperiments+1,1:obj.getnExperiments) = obj.getBasisFct{1}(para,obj.getInputData)';
%         C(obj.getnExperiments+2,1:obj.getnExperiments) = obj.getBasisFctDerivative{comb(iPara,1)}(para,obj.getInputData)';
% 
%         % Calculate the objective function
%         Cinv = inv(C);
%         saveCoeff(end-obj.BasisFctParameters+1) =(Cinv(obj.getnExperiments+1,:)*[obj.getOutputData;zeros(2,1)]);
%     end
% end



%% 
%     % Separate the different 
%     C = zeros(obj.getnExperiments+1+1);
%     C(1:obj.getnExperiments,1:obj.getnExperiments) = ...
%                         obj.CovariogramMatrix(1:obj.getnExperiments,1:obj.getnExperiments);
%     for iDerivative=1:obj.nBasisFctDerivative
%         % Basis Function
%         C(1:obj.getnExperiments,obj.getnExperiments+1) = obj.getBasisFct{1}(para,obj.getInputData);
%         % Derivatives
%         C(1:obj.getnExperiments,end) = obj.getBasisFctDerivative{iDerivative}(para,obj.getInputData);
%         % All together
%         C=triu(C)+triu(C)'-bsxfun(@times,eye(length(C)),diag(C));
% 
%         % Calculate the objective function
%         Cinv = inv(C);
% %         saveCoeff(iDerivative) =(Cinv(end,:)*[obj.getOutputData;zeros(2,1)])/sum(Cinv(end,:).^2);
%         saveCoeff(iDerivative) =(Cinv(end,:)*[obj.getOutputData;zeros(2,1)]);
%     end
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
