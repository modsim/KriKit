function [ ] = calcInverseCoVar(obj)
% Calculates the extended variance matrix and its inverse: 
% [Var(:) ones(size(Var,1),1);ones(1,size(Var,1)) 0]
% With Var(:) is the actual variance matrix with repect to the distance
% between the data points. 
% 
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.
    
    switch obj.ModelChoice
        case {1,2}
            obj.CoVar = ones(obj.nExperiments,obj.nExperiments);
%             error('calcInverseCoVar was not implemented for model 1 and 2, yet');
            for varIndex1=1:obj.nExperiments
               for varIndex2=varIndex1+1:obj.nExperiments
                   % Use the distance of the input variables between the experiment 
                   % which was already calculated. (Betw. Exp varIndex1 & varIndex2)
                   % Distance is a vector of nInputVar X 1
                   
                   % Index = 
                   Index = (varIndex1)*varIndex1/2-1+(varIndex2-varIndex1);
%                    Index = (varIndex1^2-varIndex1+2*varIndex2-2)/2;
                    obj.CoVar(varIndex1,varIndex2) =  obj.CovarModel(obj.DistInput(varIndex1+varIndex2)); % Here should be calculated the which index represent which distance
                    obj.CoVar(varIndex2,varIndex1) = obj.CoVar(varIndex1,varIndex2);
               end
            end
        case 3
            obj.CoVar = ones(obj.nExperiments,obj.nExperiments);
            for varIndex1=1:obj.nExperiments
               for varIndex2=varIndex1+1:obj.nExperiments
                   % Use the distance of the input variables between the experiment 
                   % which was already calculated. (Betw. Exp varIndex1 & varIndex2)
                   % Distance is a vector of nInputVar X 1
                    obj.CoVar(varIndex1,varIndex2) =  obj.CovarModel(obj.DistInput(:,varIndex2,varIndex1)');
                    obj.CoVar(varIndex2,varIndex1) = obj.CoVar(varIndex1,varIndex2);
               end
            end
        otherwise
    end
    
    
    
    
    obj.InvCoVar=inv(obj.CoVar);
    
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
