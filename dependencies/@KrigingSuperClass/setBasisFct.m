function [] = setBasisFct(obj,type,varargin)
% defineBasisFct(type,varargin)
% This function defines the basis functions.
%
% - 'type' is a string and can contain the elements 'polynomial' or 'user',
% ... 'polynomial' means that the basis functions are of the form :
%      F = {1,x1,x2,x1*x2,x1^2,x2^2,x1*x1^2,x1^2*x1,x1^2*x1^2,....}(if only two input variables exist)
%      "varargin" is an integer indicating the degree of the polynomial.
%      e.g.defineBasisFct('polynomial',2) leads to a basis function of the
%      form '@(p,x)x(:,1).^0+x(:,1).^1+x(:,1).^2' (case: single input
%      variable)
% ... 'user' than varagin consists of two part. First entry is an intger an provides information about the number of parameters in the function. The second is a cell containing the user defined basis function. 
% basis function should only get one input matrix! 
% The determined basis function can be checked by using getBasisFct()
% 
% E.g.: defineBasisFct('user',2,'@(p,x)p(1)*x(:,1).^2+p(2)*x(:,2)'), 
% x(:,1),x(:,2) are used since then you can make later the prediction at
% several points parallel
%
% You can get:
% - BasisFct
% - BasisFctParameters
% - BasisFctCoefficients ... Kriging coefficients of the bsis functions.
%                            Call estimateBasiisFunctionCoefficients first
% - BasisFctType ... saves "type"
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

% Check Input
if ~ischar(type)
    error('Second input of defineBasisFct should be a string!')
end
if (obj.nInputVar==0)
    warning('nInputVar is 0, you may have used "setBasisFct" before setInput')
end


% Define polynomial functions
obj.MaxDegree = [];
switch type
    case obj.allowed{1}
        % Get degree of the polynomial
        obj.MaxDegree = varargin{1};
        obj.BasisFctType = type;
        if ~isnumeric(obj.MaxDegree)||obj.MaxDegree<0
            error('First Input should be an integer bigger or equal to 0!')
        end

        % Define Polynomial
        degreeMatrix=combn(0:obj.MaxDegree,obj.nInputVar);
        degreeMatrix(sum(degreeMatrix,2)>obj.MaxDegree,:)=[];
        obj.nBasisFct = size(degreeMatrix,1);
        obj.BasisFct = cell(obj.nBasisFct,1);
        for iTerm=1:size(degreeMatrix,1)
            StringFct = '@(p,x)ones(size(x,1),1)';
%             StringFct = '@(p,x)1';
            for iInput = 1:obj.nInputVar
                StringFct=strcat(StringFct,'.*x(:,',num2str(iInput),').^',num2str(degreeMatrix(iTerm,iInput)));
            end
            eval(strcat('obj.BasisFct{iTerm}=',StringFct,';'))
        end
    case obj.allowed{end}
        % User
        if length(varargin)~=2
            error('Input in case of user should be 1. number of parameters 2. function')
        end
        obj.BasisFctType = type;

        obj.nBasisFct = 1;
        obj.BasisFct = cell(obj.nBasisFct,1);
        obj.nBasisFctParameters = varargin{1}(1);
        if obj.nBasisFctParameters > obj.nExperiments
            warning('More parameters to determine (%d) than experiments are provided (%d)',obj.nBasisFctParameters,obj.nExperiments)
        end
        obj.BasisFct{1}=eval(varargin{2});
        if ~strncmp(varargin{2},'@(p,x)',6)
            error('Incorrect form of your basis function. The provided BasisFct has to be of the form f=@(p,x)...')
        end
    otherwise
        obj.allowed
        error('Incorrect choice for "type"(Second Input). You chose: %s but allowed are only the aboved mentioned ones!',type);
end

% Test Basis Function
if nargin(obj.BasisFct{1})~=2
    error('Incorrect form of your basis function. The provided BasisFct has to be of the form f=@(p,x)...')
end

obj.BasisFct{1}(ones(obj.nBasisFctParameters,1),ones(obj.nInputVar,obj.nInputVar));
if(size(obj.BasisFct{1}(ones(obj.nBasisFctParameters,1),obj.getInputData),1)~=size(obj.getOutputData,1))
    error('size of outputData and size of basis function evaluations at the inputData must be the same. If you want have for example only one constant "a" than set "@(p,x)ones(size(x,1),1)*a" ')
end

if ~isempty(obj.getCovariogramMatrix)
    try
        if length(obj.BasisFctParameters)~=obj.nBasisFctParameters
            obj.BasisFctParameters = ones(obj.nBasisFctParameters,1);
            warning('Since number of basis function parameters changed, the parameter array is set to an array consisting of ones')
        end
        obj.calcCovariogramMatrix();
    catch exception
        warning('By redefining of the basis function the covariogram has to be recalculated. Please run Kriging.calcCovariogramMatrix after defining the basis function parameters')
        rethrow(exception)
    end
else
    if isempty(obj.BasisFctParameters)
        obj.BasisFctParameters= ones(obj.nBasisFctParameters,1);
    end
    obj.calcCovariogramMatrix
end

if obj.nBasisFct>obj.nExperiments
    error('Not enough data for estimating coefficients for all basis functions')
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
