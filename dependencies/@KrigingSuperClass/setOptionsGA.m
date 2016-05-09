function [options] = setOptionsGA(obj)
% [options] = setOptionsGA()
% This functions sets the options for the genetic algorithms
% (Called in solveLeastSquareBasisFctGA or solveLeastSquareCovariogramGA)
%
% You can set: 
% - ShowDetails ... if true, process is desplayed in command window
% - UseParallel ... if true and parallel toolbox is available, solver uses
%                   parallel process power
% - PopulationSize ... size of the population in each run
% - CR ... CrossoverFraction, fraction of the population at the next
%          generation
% - TimeLimit ... Maximum time the solver runs
% - PopInitMatrix ... Decide in which range the genes of the initial
%                     population is allowed to be. PopInitMatrix is a
%                     matrix of size nParameters x 2 
% - Generations ... maximum number of iteration after which the solver
%                   should stop. By default nIterationsSolver=1e3
%
% - UseUniformDistributedPopulation: Decide if initial population should be uniform distributed over
%                                    the entire defined parameter space.
%                                    LBBasisFctParameters and
%                                    UBBasisFctParameters must be defined
%                                    here
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

    switch obj.ShowDetails
        case 0
            showDetails='off';
        case 1
            showDetails='iter';
        otherwise
            error('ShowDetails=%i was not defined',ShowDetails);
    end
    if obj.PopulationSize<=0
        obj.PopulationSize=min(max(10*obj.nInputVar,40),100);
    end
    
    try
        if obj.getUseParallel==1
            strParallel = {'always'};
        else
            strParallel = {'never'};
        end
        options = gaoptimset('CrossoverFraction',obj.CR,...
                             'Generations',obj.Generations,...
                             'TimeLimit',obj.TimeLimit,...
                             'Display',showDetails,...
                             'PopulationSize',obj.PopulationSize,...
                             'PopInitRange',obj.PopInitMatrix',...
                             'UseParallel',strParallel{1},...
                             'EliteCount',2);
    catch ex
        
        if ~isempty(obj.PopInitMatrix)
            warning('PopInitMatrix is not correct defined')
        else
            warning('gaoptimset could not be used')
        end
        warning(ex.message)
        
        options = struct('CrossoverFraction',obj.CR,'Generations',obj.Generations,...
        'TimeLimit',obj.TimeLimit,'Display',showDetails,'PopulationSize',obj.PopulationSize);
    end
    
    if obj.UseUniformDistributedPopulation==1
        if isempty(obj.getLBBasisFctParameters)||isempty(obj.getUBBasisFctParameters)
            error('"LBBasisFctParameters" and "UBBasisFctParameters" must be defined when you want to used the option "UseUniformDistributedPopulation"')
        end
        options = gaoptimset(options,'CreationFcn',@obj.gaCreationUniformFeasible);
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
