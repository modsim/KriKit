function [bestIndividual,bestValue] = ga2(obj,fctHandle,nPara,LB,UB,options)
% ga2(fctHandle,nPara,LB,UB,options)
%
% This script implements the Simple Genetic Algorithm described
% in the examples section of the GA Toolbox manual.
%
% Input: 
%   - fctHandle ... for function evaluation (ouputMatrix = fctHandle(Input);)
%   - nPara ... Number of input variables
%   - LB/UB ... Lower and Upper of input variables
%   - options ... see documentations of "setOptionsGA()"
%
% Author:     Andrew Chipperfield
% History:    23-Mar-94     file created
%
% tested under MATLAB v6 by Alex Shenfield (22-Jan-03)
%
% Slightly modified by:
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

% Number of individuals per subpopulations
NIND = options.PopulationSize;           
% maximum Number of generations
MAXGEN = options.Generations;        
% Generation gap, how many new individuals are created
GGAP = .9;           
% nPara = 20;           % Number of variables
PRECI = 20;          % Precision of binary representation
TimeLimit = options.TimeLimit; 
if strcmp(options.Display,'iter');
    showDetails=true;
else
    showDetails=false;
end

if ( sum(isinf(LB))+sum(isinf(UB)) )>0
    warning('Lower and Upper bound but be real numbers (-)inf is not allowed and are set to (-)10^6')
    LB(isinf(LB))=-10^6;
    UB(isinf(UB))=10^6;
end

if ( sum(isnan(LB))+sum(isnan(UB)) )>0
    error('Lower and Upper bound but be real numbers (-)inf and nan are not allowed')
end

if isempty(LB)
    error('No lower bounds are defined')
end

if isempty(UB)
    error('No upper bounds are defined')
end


    if obj.ShowWaitingBar==1;
        str =  horzcat('Please wait prediction runs ... already done: ',num2str(0),'%');
        hWaitingBar = waitbar(0,str);
    end
% Build field descriptor
   try
       FieldD = [rep(PRECI,[1, nPara]); LB; UB;...
                  rep([1; 0; 1 ;1], [1, nPara])];
   catch ex
       warning('Error detected: Maybe you did not installed open source genetic algortihm tool bot (http://www.shef.ac.uk/acse/research/ecrg/gat.html)')
       rethrow ex
   end

% Initialise population
   Chrom = crtbp(NIND, nPara*PRECI);

% Reset counters
   Best = NaN*ones(MAXGEN,1);	% best in current population
   gen = 0;			% generational counter

% Evaluate initial population
   ObjV = singleApply(bs2rv(Chrom,FieldD));

% Track best individual and display convergence
   Best(gen+1) = min(ObjV);
   tElapsed = 0;

tic;   

if showDetails
    fprintf('Generation\t|\tEvaluations\t|\tBest\t|\tMean\n')
end

% Generational loop
    
    
   while gen < MAXGEN&&tElapsed < TimeLimit

%     fprintf('Current Generation: %i\n',gen);
    if obj.ShowWaitingBar==1;
        str =  horzcat('Please wait estimation runs ... already done: ',num2str(gen/MAXGEN*100,'%3.2f'),'%');
        waitbar(gen/MAXGEN,hWaitingBar,str)
    end
    % Assign fitness-value to entire population
       FitnV = ranking(ObjV);

    % Select individuals for breeding
       SelCh = select('sus', Chrom, FitnV, GGAP);
        
    % Recombine selected individuals (crossover)
       SelCh = recombin('xovsp',SelCh,0.7);

    % Perform mutation on offspring
       SelCh = mut(SelCh);

    % Evaluate offspring, call objective function
       ObjVSel = singleApply(bs2rv(SelCh,FieldD));

    % Reinsert offspring into current population
       [Chrom ObjV]=reins(Chrom,SelCh,1,1,ObjV,ObjVSel);

    % Increment generational counter
       gen = gen+1;
       tElapsed = toc;

    % Update display and record current best individual
       Best(gen+1) = min(ObjV);
       
%        plot(log10(Best),'ro'); xlabel('generation'); ylabel('log10(f(x))');
%        text(0.5,0.95,['Best = ', num2str(Best(gen+1))],'Units','normalized');
%        drawnow;

    if showDetails
        fprintf('%i|%17e|%16e|%15e|\n',gen,gen*NIND,min(ObjV),mean(ObjV))
    end
    
   end 
   
% End of GA
bestIndividual = bs2rv(Chrom(ObjV==min(ObjV),:),FieldD);
bestValue = min(ObjV);

if obj.ShowWaitingBar==1;
    close(hWaitingBar)
end
%% Nested Function
    function [ouputMatrix] = singleApply(Input)
        ouputMatrix = zeros(size(Input,1),1);
        for iIndividuals = 1:size(Input,1)
            ouputMatrix(iIndividuals) = fctHandle(Input(iIndividuals,:));
        end
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
