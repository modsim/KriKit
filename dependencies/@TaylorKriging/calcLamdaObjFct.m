function [out]=calcLamdaObjFct(obj,parameters)
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

    para = obj.getBasisFctParameters;
    k=obj.getnBasisFctParameters;
    n=obj.getnExperiments;
    
%     para = [5,2];
%     k=2;
%     n=25;
    allComb = 1+k^2+k*(k-1)/2;

    % Lambdas for all combinations of derivatives and originial basis function
    lambdaVec = parameters(1:allComb*n);
    mu1Vec = parameters(allComb*n+1:allComb*n+allComb*k);
    mu2Vec = parameters(allComb*n+allComb*k+1:allComb*n+allComb*k+k*k);

    input = obj.getInputData;
    output = obj.getOutputData;
    basis = obj.getBasisFct;
    derivatives = obj.getBasisFctDerivative;
    C=obj.getCovariogramMatrix;
    C = C(1:n,1:n);

    lambdaMatrix = reshape(lambdaVec,n,allComb);
    mu1Matrix = reshape(mu1Vec,allComb,k);
    mu2Matrix = reshape(mu2Vec,k,k);

    F = zeros(n,allComb);
    F(:,1) = basis{1}(para,input);
    for iF = 2:allComb
        F(:,iF) = derivatives{iF-1}(para,input);
    end

    ZSquare = output*output';
    % P1 
    delta1 = zeros(allComb,1);
    delta1(1+1) = 1;

    s1=C*lambdaMatrix(:,1+1) + F*mu1Matrix(:,1)-ZSquare*lambdaMatrix(:,2:k+1)*mu2Matrix(:,1);
%     s2=F'*lambdaMatrix(:,1+1)-delta1;
%     s3=2*lambdaMatrix(:,1)'*ZSquare*lambdaMatrix(:,1+k+1) - lambdaMatrix(:,1+1)'*ZSquare*lambdaMatrix(:,1+1);
%     s4=lambdaMatrix(:,1)'*ZSquare*lambdaMatrix(:,1+k+k+1) - lambdaMatrix(:,1+1)'*ZSquare*lambdaMatrix(:,1+2);

    % P1 
    delta2 = zeros(allComb,1);
    delta2(1+2) = 1;

    s5=C*lambdaMatrix(:,1+2) + F*mu1Matrix(:,2)-ZSquare*lambdaMatrix(:,2:k+1)*mu2Matrix(:,2);
%     s6=F'*lambdaMatrix(:,1+2)-delta2;
%     s7=2*lambdaMatrix(:,1)'*ZSquare*lambdaMatrix(:,1+k+2) - lambdaMatrix(:,1+2)'*ZSquare*lambdaMatrix(:,1+2);
    % s8=lambdaMatrix(:,1)'*ZSquare*lambdaMatrix(:,1+k+k+1) - lambdaMatrix(:,1+1)'*ZSquare*lambdaMatrix(:,1+2);

    s = [s1;s5];
    out=sum(s.^2);
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
