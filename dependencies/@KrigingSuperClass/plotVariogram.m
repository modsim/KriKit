function [] = plotVariogram(obj)
% This function creates a plot of the distances between the input variables
% versus variance. THe actual plot depends on your model choice and the
% number of input variable you use. 
% Case CovariogramModelChoice= 1|2 ... The plot is a 2-D plot
% Case CovariogramModelChoice= 3 & nInputVar=2 ... Three plot are generated
%   the first plot is a 3-D Plot with the absolute distances with respect to 
%   the first and second input variable on the x-axis and the y-axis,
%   respectively. On the z-axis the variance is plotted. The other two plot
%   are similar are 2-D plot which plot the variogram model setting the
%   distance with repect to the first or second input variable equal to zero.
% Case CovariogramModelChoice= 3 & nInputVar>2 Only the 2-D plots are
%   generated.
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

%     switch obj.CovariogramModelChoice
%         case {1,2,5}
    if size(obj.DistInput,2)==1
            % Plot the sample points
            figure()
            plot(obj.DistInput,obj.VarianceEstimation,'LineStyle','none','Marker','*','MarkerEdgeColor',[0.7,0.7,0.7]);
            hold on
            xData = linspace(0,max(obj.DistInput));
        %     plot(xData,obj.theta-obj.CovarModel(xData),'r','LineWidth',5);
            yData = obj.CovarModel(0,1)-obj.CovarModel(xData',0);
            plot(xData,yData,'k','LineWidth',2);
            xlabel('Distance','FontSize',20);
            ylabel('Variance','FontSize',20);
            axis([-0.1 max(obj.DistInput) min(0,min(obj.VarianceEstimation)) max(max(max(obj.VarianceEstimation)),max(yData))])
    else
        switch obj.nInputVar
            case 1
                figure()
                    % Plot measured variance
                plot(obj.DistInput,obj.VarianceEstimation,'LineStyle','none','Marker','*','MarkerEdgeColor',[0.7,0.7,0.7]);
%                     plot(obj.DistInput,obj.VarianceEstimation,'*');
                hold on
                    % plot Variogram model
                xData = linspace(0,max(obj.DistInput))';
                yData = obj.CovarModel(0,1)-obj.CovarModel(xData',0);
                plot(xData,yData,'k','LineWidth',2);
                xlabel('Distance','FontSize',20);
                ylabel('Variance','FontSize',20);
                axis([-0.1 max(obj.DistInput) min(0,min(obj.VarianceEstimation)) max(max(max(obj.VarianceEstimation)),max(yData))])
            case 2
                figure()
                hold on
                % Plot measured variance
                plot3(obj.DistInput(:,1),obj.DistInput(:,2),obj.VarianceEstimation,'*');
                    % plot Variogram model
                nGrid = 1e2;
%                     nGrid = 5;
                [x,y] = ndgrid(linspace(0,max(obj.DistInput(:,1)),nGrid),linspace(0,max(obj.DistInput(:,2)),nGrid) );
                xData = linspace(0,max(obj.DistInput(:,1)),nGrid);
                yData = linspace(0,max(obj.DistInput(:,2)),nGrid);
                zData = obj.CovarModel(zeros(1,obj.nInputVar),1)-obj.CovarModel([x(:),y(:)],0);
%                     zData = obj.CovarModel(0,1)-obj.CovarModel([y(:),x(:)],1);
                % Transpose the reshaped matrix since in the matrix the
                % rows should respresent the effect of the first input
                % variable and the column the effect of the second
                % input variable
                mesh( xData,yData,reshape(zData,sqrt(size(zData,1)),sqrt(size(zData,1)))' )
                xlabel('Distance in of input variable 1','FontSize',20)
                ylabel('Distance in of input variable 2','FontSize',20)
                zlabel('Variance','FontSize',20)

                X = linspaceNDim(zeros(1,obj.nInputVar),max(obj.DistInput,[],1),nGrid)';
                for iI = 1:obj.nInputVar
                    figure()
                    plot(obj.DistInput(:,iI),obj.VarianceEstimation,'LineStyle','none','Marker','*','MarkerEdgeColor',[0.7,0.7,0.7]);
                    hold on
                    xData = X(:,iI);
%                         distance = bsxfun(@times,ones(nGrid,obj.nInputVar),mean(obj.DistInput));
                    distance = zeros(nGrid,obj.nInputVar);
                    distance(:,iI) = X(:,iI);
                    yData = obj.CovarModel(zeros(1,obj.nInputVar),1)-obj.CovarModel(distance,0);
                    plot(xData,yData,'k','LineWidth',2);
                    xlabel('Distance','FontSize',20);
                    ylabel('Variance','FontSize',20);
                    title(strcat('Input Variable ',num2str(iI)),'FontSize',20)
                    axis([-0.1 max(xData) min(0,min(obj.VarianceEstimation)) max(max(max(obj.VarianceEstimation)),max(yData))])
                end
            otherwise
                warning('Plot Variogram for each Input variable. However, do not interprete to much in these plot since interactions between the input variable are not considered')
                nGrid = 1e1;
                X = linspaceNDim(zeros(1,obj.nInputVar),max(obj.DistInput,[],1),nGrid)';
                for iI = 1:obj.nInputVar
                    figure()
                    plot(obj.DistInput(:,iI),obj.VarianceEstimation,'LineStyle','none','Marker','*','MarkerEdgeColor',[0.7,0.7,0.7]);
                    hold on
                    xData = X(:,iI);
%                         distance = bsxfun(@times,ones(nGrid,obj.nInputVar),mean(obj.DistInput));
                    distance = zeros(nGrid,obj.nInputVar);
                    distance(:,iI) = X(:,iI);
                    yData = obj.CovarModel(zeros(1,obj.nInputVar),1)-obj.CovarModel(distance,0);
                    plot(xData,yData,'k','LineWidth',2);
                    xlabel('Distance','FontSize',20);
                    ylabel('Variance','FontSize',20);
                    title(strcat('Input Variable ',num2str(iI)),'FontSize',20)
                    axis([-0.1 max(xData) min(0,min(obj.VarianceEstimation)) max(max(max(obj.VarianceEstimation)),max(yData))])
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
