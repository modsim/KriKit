function x=marginalEI(a,b,m,s)
xNorm = bsxfun(@rdivide,bsxfun(@minus,b,m),s);
improveAbs = bsxfun(@minus,a,m);
if size(improveAbs,1)==1
    x= bsxfun(@times,gausspdf(xNorm),s) + bsxfun(@times,gausscdf(xNorm),improveAbs);
else
    x= bsxfun(@times,gausspdf(xNorm),s) + gausscdf(xNorm).*improveAbs;
end
