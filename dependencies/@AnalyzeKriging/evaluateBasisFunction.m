function []=evaluateBasisFunction(obj,varargin)
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.
    KrigingObjectIndex = varargin{1};
    Input = varargin{2};
    dimensionInterpolation = varargin{3};
    
    % Calculate linear parameters automatically applied by Kriging
    BasisFct = obj.KrigingObjects{KrigingObjectIndex}.getBasisFct;
    obj.KrigingObjects{KrigingObjectIndex}.estimateBasisFctCoefficients;
    OutputMatrix = zeros(size(Input,1),2);
    basisCoeff = obj.KrigingObjects{KrigingObjectIndex}.getBasisFctCoefficients;
    if obj.KrigingObjects{KrigingObjectIndex}.getUseSimpleKriging
        basisCoeff(:)=1;
    end

    % Actual evaluation
    for iBasisFctNested = 1 : size(BasisFct,1)
        OutputMatrix = ...
            OutputMatrix+ ... 
            basisCoeff(iBasisFctNested)*...
            [BasisFct{iBasisFctNested}(obj.KrigingObjects{KrigingObjectIndex}.getBasisFctParameters,Input),...
            zeros(size(Input,1),1)];
    end
    
    switch dimensionInterpolation
        case 2
            obj.KrigingPrediction_Interpolation2D{KrigingObjectIndex,1} = OutputMatrix;
        case 3
            obj.KrigingPrediction_Interpolation3D{KrigingObjectIndex,1} = OutputMatrix;
        otherwise
            error('Diminsion is unusual')
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
