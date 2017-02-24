function [varargout] = sampleFromDistribution(obj,varargin)
    muVec = varargin{1};
    stdVec = varargin{2};
    
    distA = randn(1e5,1);
    distB = randn(1e5,1);
    
    distA = distA*stdVec(1) + muVec(1);
    distB = distB*stdVec(2) + muVec(2);
    
    varargout{1} = distA;
    varargout{2} = distB;
end

