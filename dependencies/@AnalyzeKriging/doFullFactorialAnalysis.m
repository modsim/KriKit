function [] = doFullFactorialAnalysis(obj,KrigingObjectIndex)
% [] = doFullFactorialAnalysis(obj,KrigingObjectIndex)
% This function create a full factorial + center point design inside of the Kriging
% interpolation and applies afterwards an ANOVA. 
% 
% Input:
% - KrigingObjectIndex ... index of Kriging object of interest
% 
% Output: -
% 
% You can set:
% - ShowDetails ... a table is display which contains the variable in the of
% - LB/UBInputVarInterpolation ... determines the range of input variable
%   for the DoE Analysis
% 
% You can get: -
%
% For more information about sets and get see documentation of
% "doANOVAandPrediction".
% 
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.
%% Initialization
nInputVar = obj.KrigingObjects{KrigingObjectIndex}.getnInputVar;
inputData = obj.KrigingObjects{KrigingObjectIndex}.getInputData;
nData = 2^nInputVar;
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

factorLevels = [-1,1];

% Get Grid 
str={'['};
str2={'ndgrid('};

% Determine Area Of Interest
defineBoundOfInputVar(obj,KrigingObjectIndex);
lowerBoundsInputVariables = obj.getLBInputVarInterpolation;
lowerBoundsInputVariables = lowerBoundsInputVariables{KrigingObjectIndex};
upperBoundsInputVariables = obj.getUBInputVarInterpolation;
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

% Convert to actual values
for iVar=1:nInputVar
%     inputMatrix(:,iVar) = (normInputMatrix(:,iVar)-(-1))*(max(inputData(:,iVar))-min(inputData(:,iVar)))/(1-(-1))+min(inputData(:,iVar));
    inputMatrix(:,iVar) = (normInputMatrix(:,iVar)-(-1))*(upperBoundsInputVariables(iVar)-lowerBoundsInputVariables(iVar))/(1-(-1))+lowerBoundsInputVariables(iVar);
end

%% Do the Main Work

% Initialization
nCoeff = 1;
for iVar = 1:nInputVar
    nCoeff = nCoeff + nchoosek(nInputVar,iVar);
end
df = nData-nCoeff;
X = ones(nData,nCoeff);
indexColumnX = 2;

% Calculate Systen Matrix
for iVar = 1:nInputVar
% Initialization
combinations = nchoosek(1:nInputVar,iVar);

    for iComb = 1:size(combinations,1)
        X(:,indexColumnX) = prod(normInputMatrix(:,combinations(iComb,:)),2);
        indexColumnX = indexColumnX + 1;
    end

end
% Do Prediciton and do ANOVA based on this
doANOVAandPrediction(obj,KrigingObjectIndex,X,inputMatrix,df)

%% Output
if obj.KrigingObjects{KrigingObjectIndex}.getShowDetails==1
    fprintf('Legend: \n')
    for iVar = 1:nInputVar
        fprintf('%s := %s\n',char(iVar-1+'A'),strInputNames{iVar})
    end
    fprintf('Analysis Result: \n')
%     fprintf('%-30.18s|%17.17s|%16.16s|%15.15s|%15.15s|\n','Para Combination','Mean(Coefficient)','Std(Coefficient)','tValue','P(>|t|)')
    if isempty(obj.ANOVAPvalue)
        fprintf('%-30.18s|%17.17s|%16.16s|%15.15s|\n','Para Combination','Mean(Coefficient)','Std(Coefficient)','tValue')
        fprintf('%-30.30s|%17e|%16e|%15e|\n','Intercept',obj.ANOVACoefficients(1),obj.ANOVAStdOfCoefficients(1),obj.ANOVATvalue(1));
    else
        fprintf('%-30.18s|%17.17s|%16.16s|%15.15s|%15.15s|\n','Para Combination','Mean(Coefficient)','Std(Coefficient)','tValue','P(>|t|)')
        fprintf('%-30.30s|%17e|%16e|%15e|%15e|\n','Intercept',obj.ANOVACoefficients(1),obj.ANOVAStdOfCoefficients(1),obj.ANOVATvalue(1),obj.ANOVAPvalue(1));
    end
    indexColumnX = 2;
    % - Show Combinations
    for iVar = 1:nInputVar
    % Initialization
    combinations = nchoosek(1:nInputVar,iVar);

        for iComb = 1:size(combinations,1)
            % obj.ANOVAPvalue should be empty if no statistical toolbox is
            % installed
            if isempty(obj.ANOVAPvalue)
            fprintf('%-30.30s|%17e|%16e|%15e|\n',char(combinations(iComb,:)-1+'A'),...
                                                  obj.ANOVACoefficients(indexColumnX),...
                                                  obj.ANOVAStdOfCoefficients(indexColumnX),...
                                                  obj.ANOVATvalue(indexColumnX))
            else
            fprintf('%-30.30s|%17e|%16e|%15e|%15e|\n',char(combinations(iComb,:)-1+'A'),...
                                                  obj.ANOVACoefficients(indexColumnX),...
                                                  obj.ANOVAStdOfCoefficients(indexColumnX),...
                                                  obj.ANOVATvalue(indexColumnX),...
                                                  obj.ANOVAPvalue(indexColumnX))
            end
            indexColumnX = indexColumnX + 1;
        end

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
