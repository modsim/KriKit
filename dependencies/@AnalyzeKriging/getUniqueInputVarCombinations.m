function [uniqueRows,countUniqueRows]=getUniqueInputVarCombinations(obj,varargin)
% [uniqueRows,countUniqueRows]=getUniqueInputVarCombinations(obj,KrigingObjectIndex,indicesNotInvestigating)
%
% This function identifies input values where the highest number of sample
% points exist. 
% 
% For Input = [2,1;2,0;1,3;2,2];
% [uniqueRows,countUniqueRows]=getUniqueInputVarCombinations(1,2) will give
% out uniqueRows = [1,2]; with countUniqueRows=[1,3]; Meaning that for
% input variable 1(=not 2). There exist two different values 1 and 2 (see
% uniqueRows). Where input variable 1 take the value 1 only once but it
% takes the value 2 three times (see countUniqueRows)
% 
% Input:
% - KrigingObjectIndex ... index of the Kriging objects of interst
% - indicesNotInvestigating  ... array of indices of input variables which
%                                are continously varied and ploted on the
%                                x(and y-) axes. These input variables are
%                                not considered by identifying parameter
%                                combinations
%
% Output:
% uniqueRows ... unique combinations of input variables which are
%                non-continously varied (opposite of
%                indicesNotInvestigating)
% countUniqueRows ... Number of data points with same values as in
%                     uniqueRows 
%
% You can set: -
%
% You can get: -
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

% Initialization
KrigingObjectIndex = varargin{1};
indicesNotInvestigating = varargin{2};

InputData = obj.KrigingObjects{KrigingObjectIndex}.getInputData;
indicesInvestigating = setdiff(1:size(InputData,2),indicesNotInvestigating);

% Calculates how much parameter combination exists in the input data
% Find unique individual parameters
if isempty(indicesInvestigating)
    uniqueRows = unique(InputData,'rows');
    uniqueRows = uniqueRows(1,:);
else
    uniqueRows=unique(InputData(:,indicesInvestigating),'rows');
end


% Identify the ambundance
countUniqueRows = zeros(size(uniqueRows,1),1);
for iRow = 1:size(uniqueRows,1)
    countVector = true(size(InputData,1),1);
    for iPara = 1:length(indicesInvestigating)
        countVector = countVector&InputData(:,indicesInvestigating(iPara))==uniqueRows(iRow,iPara);
    end
    countUniqueRows(iRow) = sum(countVector);
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
