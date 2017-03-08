function [probabilityDensity]=MCMCDistributionFctPredictionSlice(obj,varargin)
% [probabilityDensity]=MCMCDistributionFctPredictionSlice(estimationpoints) 
% 
% Calculates the log(expectedImprovement) of the Kriging process at the
% estimationpoints (nXnInputVar)
%
% Input : 
% - estimationpoints ... values of input variables where the probability
%                        density function (expcted improvement) shall be
%                        evaluated (nPointsXnInputvar)
% 
% You can set:
%  - InequalityConstraintHandle ... representing a function handle to a 
%                                   function that takes a maxtrix of the
%                                   size nPointsXnInputVar as input
%                                   and gives out a nPointsX1 vector
%                                   containing true if inequality
%                                   constraints hold
% - LB(UB)InputVarInterpolation ... lower(upper) bounds of the input
%                                   variables
%
% - MinMax                      ... Decide if the output should be
%                                   maximized (+1) or minimized (-1)
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

estimationPoints = varargin{1};
nParameterSets = size(estimationPoints,1);
if ~isempty(obj.InequalityConstraintHandle)
    boolIsValid = obj.InequalityConstraintHandle(estimationPoints);
    boolIsValid = boolIsValid&~bsxfun(@lt,estimationPoints,obj.LBInputVarInterpolation{obj.ObjectiveIndicesUsedByCalcNewSamplesViaMCMC(1)});
    boolIsValid = boolIsValid&~bsxfun(@gt,estimationPoints,obj.UBInputVarInterpolation{obj.ObjectiveIndicesUsedByCalcNewSamplesViaMCMC(1)});
else
    boolIsValid = true(nParameterSets,1);
    boolIsValid = boolIsValid&~any(bsxfun(@lt,estimationPoints,obj.LBInputVarInterpolation{obj.ObjectiveIndicesUsedByCalcNewSamplesViaMCMC(1)}));
    boolIsValid = boolIsValid&~any(bsxfun(@gt,estimationPoints,obj.UBInputVarInterpolation{obj.ObjectiveIndicesUsedByCalcNewSamplesViaMCMC(1)}));
end
probabilityDensity = log(0)*ones(nParameterSets,1);

minProblemCoeff = obj.getMinMax(obj.ObjectiveIndicesUsedByCalcNewSamplesViaMCMC);
if any(boolIsValid)
    
    predOutput = obj.KrigingObjects{obj.ObjectiveIndicesUsedByCalcNewSamplesViaMCMC}.prediction(estimationPoints(boolIsValid,:));
    % Consider only the expectation value;
    predOutput = predOutput(:,1);
    % Maximization problem is converted into a minimization problem by
    % coordinate transformation.
    predOutput = minProblemCoeff*predOutput(:,1) - min(minProblemCoeff*obj.KrigingObjects{obj.ObjectiveIndicesUsedByCalcNewSamplesViaMCMC}.getOutputData);
    probabilityDensity(boolIsValid)=log(predOutput);
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
