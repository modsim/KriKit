function []=plotResultsOfOrdinaryRegressionKriging(obj)
% []=plotResultsOfOrdinaryRegressionKriging(obj)
% This function plots the results of the ordinary regression kriging when
% the number of input vairable is not bigger then 2
%
% You can set:
%
% You can get:
%
% Copyright 2014-2016: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

switch obj.nInputVar
    case 1 
        plotResultsOfOrdinaryRegressionKriging1D
    case 2
        plotResultsOfOrdinaryRegressionKriging2D
    otherwise
        error('"plotResultsOfOrdinaryRegressionKriging" was not implemented to handle dimension of input variables = %i',obj.nInputVar)
end

%% ------------------------------------------------------------------------
    function [] = plotResultsOfOrdinaryRegressionKriging1D()
        % ##### Initialization ######
        figure()
        hold on
        % ##### Plot Data #####
        plot(obj.InputData,obj.OutputData,'ko','MarkerFaceColor',[255/255 102/255 0/255]);
        % ##### Plot Interpolation #####
        plot(obj.InputData,obj.EstimatedOutputDataViaOrdinaryKriging(:,1),'x')
        % ##### Plot Regression #####
        inputData = linspace(min(obj.InputData),max(obj.InputData))';
        plot(inputData,obj.BasisFct{1}(obj.BasisFctParameters,inputData),'-');
        % ##### Labels #####
        legend('Provided Data','Kriging Estimation','Regression')
        xlabel('Input Data','FontSize',20)
        ylabel('Output Data','FontSize',20)
    end
%% ------------------------------------------------------------------------
    function [] = plotResultsOfOrdinaryRegressionKriging2D()
        % ##### Initialization ######
        figure()
        hold on
        % ##### Plot Data #####
        plot3(obj.InputData(:,1),obj.InputData(:,2),obj.OutputData,'ko','MarkerFaceColor',[255/255 102/255 0/255]);
        % ##### Plot Interpolation #####
        plot(obj.InputData(:,1),obj.InputData(:,2),obj.EstimatedOutputDataViaOrdinaryKriging,'x')
        % ##### Plot Regression #####
        calculateBasisFctOverSpace
        mesh(inputDataFirstDimensionUnique,inputDataSecondDimensionUnique,reshape(basisFunctionValues,length(inputDataFirstDimensionUnique),length(inputDataFirstDimensionUnique))' )
        % ##### Labels #####
        legend('Provided Data','Kriging Estimation','Regression')
        xlabel('Input Data','FontSize',20)
        ylabel('Output Data','FontSize',20)
        % --------------------------- Nested Function ---------------------
        function [] = calculateBasisFctOverSpace()
            nGrid = 1e2;
            % ##### Define value of the input variable at which the basis
            % function should be calculated #####
            inputDataFirstDimensionUnique = linspace(min(obj.InputData(:,1)),max(obj.InputData(:,1)),nGrid);
            inputDataSecondDimensionUnique = linspace(min(obj.InputData(:,2)),max(obj.InputData(:,2)),nGrid);
            % ##### Create Grid ##### 
            [inputDataFirstDimension,inputDataSecondDimension] = ...
                ndgrid(inputDataFirstDimensionUnique,inputDataSecondDimensionUnique);
            inputDataFirstDimension = inputDataFirstDimension(:);
            inputDataSecondDimension = inputDataSecondDimension(:);
            % ##### Calculate the values of the basis functions at the grid
            % points #####
            basisFunctionValues = obj.BasisFct{1}(obj.BasisFctParameters,[inputDataFirstDimension,inputDataSecondDimension]);
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
