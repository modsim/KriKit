function [] = plotParetoInput(obj,varargin)
%[] = plotParetoInput(krigingObjectIndex,inputIndices,plotParetoTrajectory)
%
% Plot the values of input variables which are associated to the Pareto
% front using "determineParetoSet()". Please run "determineParetoSet()"
% befor using "plotParetoInput()"
%
% Input: 
% krigingObjectIndex ... vector which contains the indices of the
%                        objectives of interest. 
%
% inputIndices ... vector which contains the indices of the input variables
%                  which distribution shall be visualized
%
% plotParetoTrajectory ... if true, input values associated with pareto
%                          optimal points are connected via lines. Input
%                          data are sorted with increasing values in
%                          objectives, starting with the objective
%                          associated with the first entry in
%                          "krigingObjectIndex". See also "sortrows"
%
% NOTE: Multiple measurement for identical input sets lead to overlapping.
%       In this case, numeration is not continuous.
%
% You can set:
% - WriteNumbers ... pareto optimal number are numerated with increasing
%                    values in the first obective.
%
% You can get: -
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

    %% Initialization
    krigingObjectIndex = varargin{1};
    inputIndices = varargin{2};
    colorNonOptimal = [158,238,255]/255;
    colorOptimal = [45,202,176]/255;
    
    if length(varargin)>2
        drawParetoTrajectory = varargin{3};
        if ~islogical(drawParetoTrajectory)
            error('"drawParetoTrajectory" has to be logical')
        end
    end
    nInputVar = length(inputIndices);
    
    % Test Input
    if (nInputVar~=1&&nInputVar~=2&&nInputVar~=3)||~any(size(inputIndices)==1)
        error('InputIndices has to be an one-dimensional array of size 1, 2, or 3')
    end
    
    if obj.nParetoSetExperiments<=0
        warning('No Pareto Set was determined. "determineParetoSet" was called')
        determineParetoSet(obj,krigingObjectIndex)
    end
    %% Plot Pareto Input
    figure();
    hold on
    grid on
    
    % Plot the sample points
    [inputSorted]=plotDataPoints();
    
    if drawParetoTrajectory
        plotParetoTrajectory(inputSorted);
    end
    
    % Label the axis
    labelAxes();
    
    % Add legend
    warning('off','MATLAB:legend:IgnoringExtraEntries')
    if nInputVar
        legend('Optimal Sample Points','Trajectory')
    else
        legend('Non-optimal Sample Points','Optimal Sample Points','Trajectory')
    end
    warning('on','MATLAB:legend:IgnoringExtraEntries')
%% Nested Functions
% ---------------------------------------------------------------------
function [inputSorted]=plotDataPoints()
    % All sample points
    inputProvided = obj.KrigingObjects{krigingObjectIndex(1)}.getInputData;
    % Sort Pareto Front Input
    [~,idxY] = sortrows(obj.ParetoSetExperiments);
    inputSorted = obj.ParetoValuesInput(idxY,:);
    [~,unqieRowIndex]=unique(inputSorted,'rows');
    
    % Actual Plotting
    switch nInputVar
    case 1
        % Mark pareto optimal
        nPareto = length(obj.ParetoValuesInput(:,inputIndices(1)));
        plot(1:nPareto,...
             inputSorted(:,inputIndices(1)),...
             'kO','MarkerFaceColor',colorOptimal,'MarkerSize',12);
        stairs(1:nPareto,...
             inputSorted(:,inputIndices(1)),...
             'k','LineWidth',2);
        set(gca,'XTick',1:nPareto)
    case 2
        if obj.ShowData==1
            % All data
            plot(inputProvided(:,inputIndices(1)),...
                 inputProvided(:,inputIndices(2)),...
                 'kO','MarkerFaceColor',colorNonOptimal,'MarkerSize',12);
        end

        % Mark pareto optimal
        plot(obj.ParetoValuesInput(:,inputIndices(1)),...
             obj.ParetoValuesInput(:,inputIndices(2)),...
             'kO','MarkerFaceColor',colorOptimal,'MarkerSize',12);
    case 3
        if obj.ShowData==1
            % All data
            plot3(inputProvided(:,inputIndices(1)),...
                  inputProvided(:,inputIndices(2)),...
                  inputProvided(:,inputIndices(3)),...
                  'kO','MarkerFaceColor',colorNonOptimal,'MarkerSize',12);
        end

        % Mark pareto optimal
        plot3(obj.ParetoValuesInput(:,inputIndices(1)),...
              obj.ParetoValuesInput(:,inputIndices(2)),...
              obj.ParetoValuesInput(:,inputIndices(3)),...
              'kO','MarkerFaceColor',colorOptimal,'MarkerSize',12);
    end
    
    % Name Sample Points if wanted
    if obj.WriteNumbers==1
        for iParetoPoint = 1:obj.nParetoSetExperiments
            if any(unqieRowIndex==iParetoPoint)
                % Define Order (From Value with low sHRR to points with high sHRR)
                switch nInputVar
                case 2
                    % Define Order (From Value with low sHRR to points with high sHRR)
                    text(inputSorted(iParetoPoint,inputIndices(1)),...
                        inputSorted(iParetoPoint,inputIndices(2)),...
                        num2str(iParetoPoint),...
                        'HorizontalAlignment','center','FontWeight','bold');
                case 3

                    text(inputSorted(iParetoPoint,inputIndices(1)),...
                         inputSorted(iParetoPoint,inputIndices(2)),...
                         inputSorted(iParetoPoint,inputIndices(3)),...
                         num2str(iParetoPoint),...
                         'HorizontalAlignment','center','FontWeight','bold');
                end
            end
        end
    end
    
end
% ---------------------------------------------------------------------
function []=plotParetoTrajectory(inputSorted)
    
    switch nInputVar
    case 2
        plot(inputSorted(:,inputIndices(1)),inputSorted(:,inputIndices(2)),'LineWidth',2);
    case 3
        plot3(inputSorted(:,inputIndices(1)),inputSorted(:,inputIndices(2)),inputSorted(:,inputIndices(3)),'LineWidth',2);
    end
    
end
% ----------------------------------------------------------------------
    function []=labelAxes()
        if isempty(obj.InputVarNames{krigingObjectIndex(1)})
            switch nInputVar
                case 1
                    xlabel('Point Number','FontSize',20)
                    xlabel(horzcat('Input Variable ',num2str(inputIndices(1))),'FontSize',20);
                case {2,3}
                    xlabel(horzcat('Input Variable ',num2str(inputIndices(1))),'FontSize',20);
                    ylabel(horzcat('Input Variable ',num2str(inputIndices(2))),'FontSize',20);
                    if nInputVar>2
                        zlabel(horzcat('Input Variable ',num2str(inputIndices(3))),'FontSize',20);
                    end
            end
        else
            
            switch nInputVar
                case 1
                    xlabel('Point Number','FontSize',20)
                    inputNames = obj.getInputVarNames(krigingObjectIndex(1));
                    ylabel(inputNames(inputIndices(1)),'FontSize',20)
                case {2,3}
                    xlabel(obj.InputVarNames{krigingObjectIndex(1)}(inputIndices(1)),'FontSize',20);
                    ylabel(obj.InputVarNames{krigingObjectIndex(1)}(inputIndices(2)),'FontSize',20);
                    if nInputVar>2
                        zlabel(obj.InputVarNames{krigingObjectIndex(1)}(inputIndices(3)),'FontSize',20);
                    end
            end
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
