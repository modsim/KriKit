function [X] = createNDGRID(LB,UB,nbinsEachDimension)
% [X] = createNDGRID(LB,UB,nbinsEachDimension)

    % MakeGrid
    nVariables = length(UB);
    XLin = linspaceNDim(LB,UB,nbinsEachDimension)';
    
    if nVariables>1
        X = [];
        for iX = 1:nVariables-1
            if iX==1
                [x1,xi] = ndgrid(XLin(:,iX),XLin(:,iX+1));
                X = [x1(:),xi(:)];
            else
                [x1,xi] = ndgrid(X(:,1),XLin(:,iX+1));
                X = [repmat(X,nbinsEachDimension,1),xi(:)];
            end

        end
    else
        X = XLin';
    end
end

