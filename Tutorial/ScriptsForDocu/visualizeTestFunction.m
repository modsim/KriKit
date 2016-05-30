function [] = visualizeTestFunction(varargin)
    
    % Generate colormap used by KriKit
    figure
    wMap = winter();
    wMap = wMap(end:-1:1,:);
    close(gcf)

    % Generate Data on 100X100 grid in the range [0,10]
    nLevelsEachDimension = 1e2;
    inputProto = createNDGRID(zeros(2,1),ones(2,1)*10,nLevelsEachDimension);
    
    % Show Results for setting one input variable fix at one
    input = ones(nLevelsEachDimension^2,3);
    
    % Create Plots
%     subplot(2,2,1)
    for iFix = sort(1:3,'descend')
%         position = sort(1:3,'descend');
        indicesNotFixed = setdiff(1:3,iFix);
        
        % Create Final Grid
        input(:,indicesNotFixed) = inputProto;
        input(:,iFix) = 2;
        
        % Calculate Output
        output = tutorialFunction(input);
        
        % Plot Result
%         subplot(2,2,position(iFix))
        figure
        mesh(unique(input(:,indicesNotFixed(1))),...
            unique(input(:,indicesNotFixed(2))),...
            reshape(output,nLevelsEachDimension,nLevelsEachDimension)')
        colormap(gcf,wMap)
        xlabel(horzcat('Input ',num2str(indicesNotFixed(1))))
        ylabel(horzcat('Input ',num2str(indicesNotFixed(2))))
        zlabel('Output')
        grid on
        alpha 0.75
        
        % Set format appropriate for documentation
        set(gcf,'Position', [100, 100, 330, 330/4*3]);
        set(gca,'FontSize',10)
        
    end
    
%     set(gcf,'Position', [100, 100, 1000, 750]);
end
