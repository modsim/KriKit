function [extX,extY]=manhattan(X,Y,sortFlag)
    
%% Sort if not sorted already
    if(sortFlag)
        [sX,I]= sort(X);
        sY = Y(I);
        extX = zeros(2*length(X)-1,1);
        extY = zeros(2*length(Y)-1,1);
    else
        sX = X;
        sY = Y;
    end
    
%% Calculate the required vector
    for iX=1:length(X)
        
        % X
        extX(iX+iX-1) = sX(iX);
        if(iX<length(X))
            extX(iX+iX) = sX(iX);
        end
        
        % Y
        extY(iX+iX-1) = sY(iX);
        if(iX<length(X))
            extY(iX+iX) = sY(iX+1);
        end
        
    end
    
end