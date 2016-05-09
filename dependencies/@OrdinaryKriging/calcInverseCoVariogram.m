function [ ] = calcInverseCoVariogram(obj)
% Calculates the extended variance matrix and its inverse: 
% [Var(:) ones(size(Var,1),1);ones(1,size(Var,1)) 0]
% With Var(:) is the actual variance matrix with repect to the distance
% between the data points. 
% 
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.
    obj.CoVariogram = ones(obj.nExperiments+1,obj.nExperiments+1)*obj.sigma;
    for(varIndex1=1:obj.nExperiments)
       for(varIndex2=varIndex1+1:obj.nExperiments)
           % Use euclidean
%            fprintf('v1: %i, v2: %i, v1-v2+1:%i v1-v2+1==1: %i\n',varIndex1,varIndex2,varIndex1-varIndex2+1,(varIndex1-varIndex2+1)==1);
            obj.CoVariogram(varIndex1,varIndex2) = obj.CovarModel(norm(obj.InputData(varIndex1,:)-obj.InputData(varIndex2,:),2),varIndex1-varIndex2+1);
            obj.CoVariogram(varIndex2,varIndex1) = obj.CoVariogram(varIndex1,varIndex2);
       end
    end
    
    % Add the sigmaError for the diagonals (To allow multiple measurements)
    obj.CoVariogram = obj.CoVariogram;
    obj.CoVariogram = obj.CoVariogram-eye(size(obj.CoVariogram,1))*(obj.sigmaError);
    obj.CoVariogram(:,end) = [ones(obj.nExperiments,1);0];
    obj.CoVariogram(end,:) = [ones(1,obj.nExperiments,1),0];
    obj.InvCoVariogram=inv(obj.CoVariogram);
    
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
