function []=determineParetoSet(obj,varargin)
% []=determineParetoSet(krigingObjIndices)
% Determine Pareto set from data base and save the result in 
%
% Input: 
% krigingObjectIndex ... vector which contains the indices of the
%                        objectives of interest. 
%
% You can set:
% - MinMax ... indicates if the optimization goal is to minimize or
%              maximize the particular objective
%
% You can get:
% - nParetoSetExperiments ... number of pareto optimal points
% - ParetoSetExperiments ... matrix containing the pareto optimal
%                            objective values.
%                            (nParetoSetExperimentsXnKrigingObjective)
% - ParetoValuesInput ... Matrix containing the input values of the
%                            associated with the pareto optimal objective
%                            values (nParetoSetExperimentsXnInputVar)
% 
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

%% Intialization
    krigingObjIndices = varargin{1};
    nObjectives = length(krigingObjIndices);
    
    % Assume that all output generated for the same input values
    nOutputs = obj.KrigingObjects{krigingObjIndices(1)}.getnExperiments;
    outputMatrix = zeros(nOutputs,nObjectives);
    inputMatrix = obj.KrigingObjects{krigingObjIndices(1)}.getInputData;
    
    
    MinMax = zeros(1,nObjectives);
    for iObj=1:nObjectives
        outputMatrix(:,iObj) = obj.KrigingObjects{krigingObjIndices(iObj)}.getOutputData;
        MinMax(iObj) = obj.getMinMax(krigingObjIndices(iObj));
    end
    
%% Actual calculation
    % Convert a Maximization Problem into a minimization problem
    [obj.ParetoSetExperiments]=determineParetoSet_Mex(bsxfun(@times,outputMatrix,-MinMax));
    obj.ParetoSetExperiments(:,:) = bsxfun(@times,obj.ParetoSetExperiments,-MinMax);
    booleanPartOfParetoSet=true(nOutputs,1);
    
    for iPoint=1:nOutputs
        for iObj=1:nObjectives
            booleanPartOfParetoSet(iPoint) = booleanPartOfParetoSet(iPoint)&any(obj.ParetoSetExperiments(:,iObj)==outputMatrix(iPoint,iObj));
        end
    end
    
%% Save Outcome
    obj.ParetoValuesInput = inputMatrix(booleanPartOfParetoSet,:);
    obj.nParetoSetExperiments = size(obj.ParetoSetExperiments,1);
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
