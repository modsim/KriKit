function [newSamplePoints] = calcNewSamplesViaMCMC(obj,varargin)
% [newSamplePoints] = calcNewSamplesViaMCMC(krigingObjectIndices,algorithm)
%
% New Samples are determined following a distribution based on the expected
% improvement of the approximated Gaussian Process Regression (Kriging)
% model. The distribution is estimated via a Monte Carlo Markov Chain
% approach using one of the following algorithms:
% 1. DRAM
% 2. Slice Sampling
% DRAM is implemented in the MCMC toolbox for Matlab. Available at
% http://helios.fmi.fi/~lainema/mcmc/ (access date 10/02/2016)
% Slice Sampling is implemented in the statistic Matlab Toolbox
%  Inequality constraints can be applied using the member variable
% 'InequalityConstraintHandle' representing a function handle to a function
% that takes a nXnInputVar matrix as input and gives out a nX1 vector
%
% Input:
% - krigingObjectIndices ... vector containing the indices of the kriging
%                            objectives of interest (1XnObjectives)
% - algorithm ... string describing the appleed MCMC approach ('DRAM' or
%                 'Slice') 
% 
% You can set:
% - nMCMCLinks ... number of links calculate for the Markov chain
% - nNewSamples ... number of samples randomly drawn from the chain
% - nCutLinks ... number of samples cutted from the beginning of the
%                 calculated Markov chain before samples are drawn
%  
%  You can get: -
%
%
% DRAM is implemented in the MCMC toolbox for Matlab. Available at
% http://helios.fmi.fi/~lainema/mcmc/ (access date 10/02/2016)
%
% For further information about setting set documentation of
% "MCMCDistributionFctDRAM()" and "MCMCDistributionFctSlice()"
% 
% 
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

% Read Input
krigingObjIndexVec = varargin{1};
obj.ObjectiveIndicesUsedByCalcNewSamplesViaMCMC = krigingObjIndexVec;

% Assume constant number of input variables and same ranges everywhere
firstKrigingIndex = krigingObjIndexVec(1);
nInputVar = obj.KrigingObjects{firstKrigingIndex}.getnInputVar;
defineBoundOfInputVar(obj,krigingObjIndexVec);

if length(varargin)>=2
    MCMCAlgorithm = varargin{2};
    switch MCMCAlgorithm 
        case 'DRAM'
        case 'Slice'
        otherwise
            error('Please choose "DRAM" or "Slice" as MCMCAlgorithm')
    end
else
    MCMCAlgorithm = 'DRAM';
end


switch MCMCAlgorithm 
    case 'DRAM'
        modelPara = {}; % Not necessary here since no model parameters exist

        % Initialize Input Variables
        inputVarProp = cell(1,nInputVar);
        for iInputVar=1:nInputVar
            % UB = obj.UBInputVarInterpolation{firstKrigingIndex}(iInputVar);
            % LB = obj.LBInputVarInterpolation{firstKrigingIndex}(iInputVar);
            % randIniSample = rand(1)*(UB-LB)+LB;
            randIniSample = rand(1)*(obj.UBInputVarInterpolation{firstKrigingIndex}(iInputVar) - ...
                                     obj.LBInputVarInterpolation{firstKrigingIndex}(iInputVar)) + ...
                                     obj.LBInputVarInterpolation{firstKrigingIndex}(iInputVar);
            inputVarProp{iInputVar} = {sprintf('x_%d',iInputVar),...
                                       randIniSample,...
                                       obj.LBInputVarInterpolation{firstKrigingIndex}(iInputVar),...
                                       obj.UBInputVarInterpolation{firstKrigingIndex}(iInputVar)};
        end

        % Define Expected improvement as distribution function
        model.ssfun     = @obj.MCMCDistributionFctDRAM;
        options.method  = 'dram';
        options.waitbar = false;
        options.verbosity = false;
        if obj.getnMCMCLinks<=0
            error('nMCMCLinks has to be bigger than 0')
        end
        options.nsimu   = obj.getnMCMCLinks;

        [~,sampleMatrix,~,sschain] = mcmcrun(model,modelPara,inputVarProp,options);
        sampleMatrixCutted = sampleMatrix(obj.nCutLinks+1:end,:);
        if obj.ConsiderOnlyMaxExpectedImprovement
            % sschain contains -2log(expectedImprovement)
            [~,idx]= sort(exp(sschain(obj.nCutLinks+1:end)/-2),'descend');
            chosenIndex = idx(1:obj.getnNewSamples);
        else
            chosenIndex = randi(size(sampleMatrixCutted,1),obj.getnNewSamples,1);
        end
        newSamplePoints = sampleMatrixCutted(chosenIndex,:);

    case 'Slice'
        % Initialize Input Variables
        randIniSample = rand(1,nInputVar).*(obj.UBInputVarInterpolation{firstKrigingIndex} - ...
                        obj.LBInputVarInterpolation{firstKrigingIndex}) + ...
                        obj.LBInputVarInterpolation{firstKrigingIndex};

        % Define MCMC Distribution function
        distributionFct = @obj.MCMCDistributionFctSlice;

        % Do Slice Sampling
        sampleMatrix = slicesample(randIniSample,...
                                   obj.getnMCMCLinks,...
                                   'logpdf',distributionFct,...
                                   'width',((obj.UBInputVarInterpolation{firstKrigingIndex}-obj.LBInputVarInterpolation{firstKrigingIndex})/5));

        % Cutting out
        sampleMatrixCutted = sampleMatrix(obj.nCutLinks+1:end,:);

        if obj.ConsiderOnlyMaxExpectedImprovement
            error('ConsiderOnlyMaxExpectedImprovement cannot be used with sliceSampling')
        else
            % Random Sampling according to approximimated distribution
            chosenIndex = randi(size(sampleMatrixCutted,1),obj.getnNewSamples,1);
        end

        newSamplePoints = sampleMatrixCutted(chosenIndex,:);
end
    
    % Only temporary
    obj.ObjectiveIndicesUsedByCalcNewSamplesViaMCMC = [];
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
