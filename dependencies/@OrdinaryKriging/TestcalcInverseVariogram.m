function [ ] = calcInverseVariogram(obj)
% Calculates the extended variance matrix and its inverse: 
% [Var(:) ones(size(Var,1),1);ones(1,size(Var,1)) 0]
% With Var(:) is the actual variance matrix with repect to the distance
% between the data points. 
% 
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.
    obj.Variogram = zeros(obj.nExperiments+1,obj.nExperiments+1);
    for(varIndex1=1:obj.nExperiments)
       for(varIndex2=varIndex1+1:obj.nExperiments)
           % Use euclidean
            delta = norm(obj.InputData(varIndex1,:)-obj.InputData(varIndex2,:),2);
            obj.Variogram(varIndex1,varIndex2) =  obj.CovarModel(0,1) - obj.CovarModel(delta,(varIndex1-varIndex2)==1);
            obj.Variogram(varIndex2,varIndex1) = obj.Variogram(varIndex1,varIndex2);
       end
    end
    
    
%     obj.Variogram = obj.Variogram-(obj.sigmaError);
%     obj.Variogram = obj.Variogram-eye(size(obj.Variogram,1))*(obj.sigmaError);
    obj.Variogram(:,end) = [ones(obj.nExperiments,1);0];
    obj.Variogram(end,:) = [ones(1,obj.nExperiments,1),0];
    obj.InvVariogram=inv(obj.Variogram);
    
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
