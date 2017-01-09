function normEI=gaussEI(mu,sigma,lb,ub)
    normVar = (ub-mu)./sigma;
    normEI = sigma.*gausspdf(normVar) + (lb-mu).*gausscdf(normVar);
end