function [X] = createNDGRID(varargin)
% [X] = createNDGRID(LB,UB,nbinsEachDimension,extraPoints)
    LB = varargin{1};
    UB = varargin{2};
    nbinsEachDimension = varargin{3};
    if length(varargin)>3
        extraPoints = varargin{4};
    else
        extraPoints = [];
    end


    % MakeGrid
    nVariables = length(UB);
    XLin = linspaceNDim(LB,UB,nbinsEachDimension)';
    XLin = ([XLin;extraPoints]);
    
    if nVariables>1
        X = [];
        for iX = 1:nVariables-1
            % It is important to chck uniqueness if extra points are
            % applied
            if iX==1
                X1 = unique(XLin(:,iX));
                X2 = unique(XLin(:,iX+1));
                
                [x1,xi] = ndgrid(X1,X2);
                X = [x1(:),xi(:)];
            else
                X1 = (X(:,1));
                X2 = unique(XLin(:,iX+1));
                
                [x1,xi] = ndgrid(X1,X2);
                X = [repmat(X,size(xi,2),1),xi(:)];
            end

        end
    else
        X = XLin';
    end
    
end

