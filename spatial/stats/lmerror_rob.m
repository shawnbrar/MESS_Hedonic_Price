function result = lmerror_rob(y,x,W);
% PURPOSE: LM error robust statistic for spatial error autocorrelation in the presence of
%  a spatially lagged dependent variable 
% ---------------------------------------------------
%  USAGE: result = lmerror_rob(y,x,W)
%  where: y = dependent variable vector
%         x = independent variables matrix
%         W = contiguity matrix (standardized)
% ---------------------------------------------------
%  RETURNS: a structure variable
%         result.meth = 'lmerror_rob'
%         result.lmerror_rob   = LM statistic
%         result.prob = marginal probability
%         result.chi1 = 6.635 (chi-squared 1 dof at 99% level)
%         result.nobs = # of observations
%         result.nvar = # of variables
% ---------------------------------------------------
% NOTE: lm > 6.635,  => small prob,
%                    => reject HO: of no spatial correlation
% ---------------------------------------------------
% See also:  lmerror, lmlag, lmlag_rob, walds, lratios, moran
% ---------------------------------------------------
% REFERENCES: Anselin, Bera, Florax, Yoon (RSUE, 1996) pages 82-83.
%             Anselin and Bera (1998) pages page 276.
%             Florax, Folmer, Rey (forthcoming RSUE, 2003)
%             (typos corected).
% ---------------------------------------------------


if nargin ~= 3
error('Wrong # of arguments to lmerror_rob');
end;
[n k] = size(x);
% do ols to get residuals
b = x\y; e = y-x*b;
epe = (e'*e)/n;

I = eye(n);
m = I-(x*inv(x'*x)*x');
g = trace((W+W')*W);
j = ((W*x*b)'*m*(W*x*b)+g*epe)/(epe*n);
num = (e'*W*e)/epe-(g*e'*W*y)/(n*j*epe);
denom = g*(1-g/(n*j));
lmerror_rob = (num*num)/denom;

prob = 1-chis_prb(lmerror_rob,1);

result.meth = 'lmerror_rob';
result.lmerror_rob = lmerror_rob;
result.prob = prob;
result.chi1 = 6.635;
result.nobs = n;
result.nvar = k;
