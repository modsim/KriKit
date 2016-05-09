function y=betabinpf(x,n,a,b)
% BETABINPF Beta Binomial probability function
% BETABINPF(x,n,a,b)

% Mean n*a/(a+b)
% Var  n*a*b(n+a+b)/(a+b)^2/(1+a+b)

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 1.2 $  $Date: 2012/09/27 11:47:33 $
y = binom(n,x).*beta(x+a,n-x+b)./beta(a,b);
