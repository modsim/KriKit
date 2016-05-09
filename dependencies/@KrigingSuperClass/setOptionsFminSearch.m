function [options] = setOptionsFminSearch(obj)
% [options] = setOptionsFminSearch()
% This functions sets the options which are used for the fmincon function
% of Matlab (Called in solveLeastSquareBasisFct and
% solveLeastSquareCovariogram)
%
% You can Set:
% - showDetails ... if true, process is desplayed in command window
% - UseParallel ... if true and parallel toolbox is available, solver uses
%                   parallel process power
% - nIterationsSolver ... maximum number of iteration after which the
%                         solver should stop. By default
%                         nIterationsSolver=1e3  
%
% You can get: -
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
   
    if obj.getUseParallel==1
            strParallel = {'always'};
        else
            strParallel = {'never'};
    end
        
    options = optimset('Display',showDetails,'UseParallel',strParallel{1},'MaxIter',obj.nIterationsSolver);
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
