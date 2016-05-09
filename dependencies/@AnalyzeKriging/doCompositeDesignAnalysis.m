function [] = doCompositeDesignAnalysis(obj,KrigingObjectIndex)
% [] = doCompositeDesignAnalysis(obj,KrigingObjectIndex)
% This function create a central composite design inside of the Kriging
% interpolation and applies afterwards an ANOVA. 
% KrigingObjectIndex ... Index of the Kriging object on which the ANOVA
% should be applied.
% You can get:
% ANOVACoefficients ... contains the polynomial coefficient calculated during the
%                       ANOVA-Analysis
% ANOVACovMatrix ... Contain the covariance of the polynomial coefficient calculated
%                    during the ANOVA-Analysis
% ANOVAStdOfCoefficients ... contains the standard deviation of the polynomial coefficient
%                            calculated during the ANOVA-Analysis
% ANOVATvalue ... contain the t-value for t-test in order to decide which of the
%                 polynomial coefficient differ significantly from zero
% ANOVAPvalue ... Contain the p-value for t-test in order to decide which of the
%                 polynomial coefficient differ significantly from zero. Ususally
%                 values lower than 0.05 represent significant coefficients
% You can set:
% ShowDetails ... a table is display which contains the variable in the of
% "You can get"
% Be aware that all get variable are overwritten everytime you run an
% ANOVA!
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

%% Initialization
nInputVar = obj.KrigingObjects{KrigingObjectIndex}.getnInputVar;
inputData = obj.KrigingObjects{KrigingObjectIndex}.getInputData;
nData = 2^nInputVar+2*nInputVar+1;
normInputMatrix = zeros(nData,nInputVar);
inputMatrix = zeros(nData,nInputVar);


%% Experimental design

if isempty(obj.getInputVarNames(KrigingObjectIndex))
    strInputNames = cell(nInputVar,1);
    for iVar = 1:nInputVar
        strInputNames{iVar}=strcat('InputVar',num2str(1));
    end
else
    strInputNames = obj.getInputVarNames(KrigingObjectIndex);
end

alphaValue = 2^(nInputVar/4);
factorLevels = [-1,1];

% Get Grid 
str={'['};
str2={'ndgrid('};

% Determine Area Of Interest
defineBoundOfInputVar(obj,KrigingObjectIndex);
lowerBoundsInputVariables = obj.getLBInputVarInterpolation;
upperBoundsInputVariables = obj.getUBInputVarInterpolation;
if ~iscell(lowerBoundsInputVariables)
    lowerBoundsInputVariables = {lowerBoundsInputVariables};
end
if ~iscell(upperBoundsInputVariables)
    upperBoundsInputVariables = {upperBoundsInputVariables};
end
lowerBoundsInputVariables = lowerBoundsInputVariables{KrigingObjectIndex};
upperBoundsInputVariables = upperBoundsInputVariables{KrigingObjectIndex};

for iVar=1:nInputVar
    str=strcat(str,'q',num2str(iVar),',');
    str2=strcat(str2,'factorLevels,');
end
str=strcat(str{1}(1:end-1),']=',str2{1}(1:end-1),');');
eval(str);
for iVar=1:nInputVar
    eval(strcat('normInputMatrix(1:',num2str(2^nInputVar),',iVar)=q',num2str(nInputVar-iVar+1),'(:);'))
end

% Star Points
for iVar=1:nInputVar
    normInputMatrix(2^nInputVar+2*(iVar-1)+1:2^nInputVar+2*(iVar-1)+2,iVar) = [alphaValue;-alphaValue];
end

% Convert to actual values
for iVar=1:nInputVar
%     inputMatrix(:,iVar) = (normInputMatrix(:,iVar)-(-1))*(max(inputData(:,iVar))-min(inputData(:,iVar)))/(1-(-1))+min(inputData(:,iVar));
    inputMatrix(:,iVar) = (normInputMatrix(:,iVar)-(-1))*(upperBoundsInputVariables(iVar)-lowerBoundsInputVariables(iVar))/(1-(-1))+lowerBoundsInputVariables(iVar);
end

%% Do the Main Work
% Initialization
combinations = nchoosek(1:nInputVar,2);
nComb = nchoosek(nInputVar,2);
nPara = 1+2*nInputVar+nComb;
df = nData-nPara;
X = zeros(nData,1+2*nInputVar+nComb);

% Calculate Systen Matrix
% - Linear Main Effect
X(:,1:1+nInputVar) = [ones(nData,1),normInputMatrix];
% - Second Order Interactions
for iComb = 1:nComb
    X(:,1+nInputVar+iComb) = prod(normInputMatrix(:,combinations(iComb,:)),2);
end
% - Quadratic Main effect
X(:,1+nInputVar+nComb+1:end) = normInputMatrix.^2;

% Do Prediciton and do ANOVA based on this
doANOVAandPrediction(obj,KrigingObjectIndex,X,inputMatrix,df)

%% Output
if obj.KrigingObjects{KrigingObjectIndex}.getShowDetails==1
    fprintf('Analysis Result: \n')
    fprintf('%-30.18s|%17.17s|%16.16s|%15.15s|%15.15s|\n','Para Combination','Mean(Coefficient)','Std(Coefficient)','tValue','P(>|t|)')
    fprintf('%-30.30s|%17e|%16e|%15e|%15e|\n','Intercept',obj.ANOVACoefficients(1),obj.ANOVAStdOfCoefficients(1),obj.ANOVATvalue(1),obj.ANOVAPvalue(1));
    % - Linear Main Effect
    for iVar = 1:nInputVar
        fprintf('%-30.30s|%17e|%16e|%15e|%15e|\n',strInputNames{iVar},obj.ANOVACoefficients(1+iVar),obj.ANOVAStdOfCoefficients(1+iVar),...
                                                  obj.ANOVATvalue(1+iVar),obj.ANOVAPvalue(1+iVar))
    end

    % - First Interaction Effect
    for iComb = 1:nComb
        fprintf('%-30.30s|%17e|%16e|%15e|%15e|\n',strcat(strInputNames{combinations(iComb,1)},'*',strInputNames{combinations(iComb,2)}),...
                                                  obj.ANOVACoefficients(1+nInputVar+iComb),obj.ANOVAStdOfCoefficients(1+nInputVar+iComb),...
                                                  obj.ANOVATvalue(1+nInputVar+iComb),obj.ANOVAPvalue(1+nInputVar+iComb))
    end

    % - Quadratic Main Effect
    for iVar = 1:nInputVar
        fprintf('%-30.30s|%17e|%16e|%15e|%15e|\n',strcat(strInputNames{iVar},'^2'),obj.ANOVACoefficients(1+nInputVar+nComb+iVar),...
                                                  obj.ANOVAStdOfCoefficients(1+nInputVar+nComb+iVar),obj.ANOVATvalue(1+nInputVar+nComb+iVar),...
                                                  obj.ANOVAPvalue(1+nInputVar+nComb+iVar))
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
