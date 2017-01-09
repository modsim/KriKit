function [varargout] = estimateBasisFctParametersFsolve(obj)
%
%
% You can set:
%  - nMultiStartPoints .... number points for applying fsolve
%  - LBBasisFctParameters/UBBasisFctParameters
% 
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.
    
    if obj.NormOutput||obj.NormInput
        error('estimateBasisFctParametersFsolve cannot be used for normalized input or/and output')
    end
    %% Initialization
    options = optimoptions('fsolve');
    options.Display='off';
%     fTestFct=obj.KrigingObj;
    if isempty(obj.getLBBasisFctParameters)||isempty(obj.getUBBasisFctParameters)
        error('LBBasisFctParameters and UBBasisFctParameters have to be defined')
    else
        Pstart = createNDGRID(obj.getLBBasisFctParameters,obj.getUBBasisFctParameters,obj.nMultiStartPoints);
    end
    
    % Allocating Memory
    Pend = zeros((obj.nMultiStartPoints)^2,obj.getnBasisFctParameters);
    qualityMeasure = ones((obj.nMultiStartPoints)^obj.getnBasisFctParameters,1);
    testVec = ones((obj.nMultiStartPoints)^obj.getnBasisFctParameters,obj.getnBasisFctParameters);
    
    %% Actual Estimation
    for iStart = 1:(obj.nMultiStartPoints)^2
        Pend(iStart,:) = fsolve(@obj.calcKrigingParaEstimationObjFctVector,Pstart(iStart,:),options);
        qualityMeasure(iStart) = sum(abs(obj.calcKrigingParaEstimationObjFctVector(Pend(iStart,:))));
        testVec(iStart,:) =          abs(obj.calcKrigingParaEstimationObjFctVector(Pend(iStart,:)))';
    end
    
    %   Consider only feasible results
    idxLBValid = all(bsxfun(@ge,Pend,obj.getLBBasisFctParameters),2);
    idxUBValid = all(bsxfun(@le,Pend,obj.getUBBasisFctParameters),2);
    idxValid = find(idxLBValid&idxUBValid);
    [~,idxMin]=min(sum(testVec(idxValid,:),2));
    minQualValid = idxValid(idxMin);
    pOpt = Pend(minQualValid,:);
    
    %% Save Results
    if ~isempty(pOpt)
        varargout{1} = pOpt;
        varargout{2} = sum((obj.calcKrigingParaEstimationObjFctVector(pOpt)).^2);
    else
        varargout{1} = zeros(1,2);
        varargout{2} = inf;
        warning('no feasible parameter set')
    end
    
    obj.setBasisFctParameters(varargout{1});
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
