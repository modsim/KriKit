function [modEI]=modExpImprByOutputConstr(obj,varargin)
% [modEI]=modExpImprByOutputConstr(originalEI,estimationMean,estimationSD) 
% 
% Input : 
% 
% 
% You can set:
% - DegreeOfOutputContraint ... adjusting weight of inequality contraints
%                               for the output. The heigher more important
%                               is are the constraints
% - InequalityConstraintOutputDistribution ... (only in case of BayesianOptimization)
%            A function handle of the form
%            [constrainsMean,constrainsSD] =  inequalityConstraintOutput(estimationMean,estimationSD)
%            With
%            estimationMean ... kriging prediction at the point of interest
%                              given in "estimationPoints" [nPointsXnObj]
%            estimationSD ... kriging prediction error at the point of
%                             interest, given in "estimationPoints"
%                             [nPointsXnObj] 
%            constrainsMean ... expected value for the contrains
%                               [nPointsXnConstrains] 
%            constrainsMean ... standard deviation for the contrains
%                               [nPointsXnConstrains] 
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

originalEI = varargin{1};
estimationMean = varargin{2};
estimationSD = varargin{3};
    
handleConstr = obj.InequalityConstraintOutputDistribution;
[contraintMean,constraintSD] = handleConstr(estimationMean,estimationSD);
probValidConstraint = Phi(-contraintMean./constraintSD);

modEI = originalEI.*(prod(probValidConstraint,2)^obj.DegreeOfOutputContraint);

%% Nested Function
function [cumProb]=Phi(x)
    cumProb = (1/2+1/2*erf(x/sqrt(2)));
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
