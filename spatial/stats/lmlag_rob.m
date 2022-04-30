function result = lmlag_rob(y,x,W);
% PURPOSE: LM lag robust statistic for a spatially dependent variable in the presence of
% spatial error autocorrelation
% ---------------------------------------------------
%  USAGE: result = lmlag(y,x,W)
%  where: y = dependent variable vector
%         x = independent variables matrix
%         W = contiguity matrix (standardized)
% ---------------------------------------------------
%  RETURNS: a structure variable
%         result.meth = 'lmlagrob'
%         result.lmlagrob   = LM statistic
%         result.prob = marginal probability
%         result.chi1 = 6.635 (chi-squared 1 dof at 99% level)
%         result.nobs = # of observations
%         result.nvar = # of variables
% ---------------------------------------------------
% NOTE: lm > 6.635,  => small prob,
%                    => reject HO: of no spatial correlation
% ---------------------------------------------------
% See also:  lmlag, lmerror, lmerror_rob, walds, lratios, moran
% ---------------------------------------------------
% REFERENCES: Anselin, Bera, Florax, Yoon (RSUE, 1996) pages 83-84.
%             Anselin and Bera (1998) page 277.
%             Florax, Folmer, Rey (forthcoming RSUE, 2003)
%             (typos corected)
% ---------------------------------------------------



if nargin ~= 3
error('Wrong # of arguments to lmlagrob');
end;
[n k] = size(x);
% do ols to get residuals
b = x\y; e = y-x*b;
epe = (e'*e)/n;

I = eye(n);
m = I-(x*inv(x'*x)*x');
g = trace((W+W')*W);
j = ((W*x*b)'*m*(W*x*b)+g*epe)/(epe*n);
lm3 = (e'*W*y)/epe-(e'*W*e)/epe;
lmlag_rob = (lm3*lm3)*(1/((n*j)-g));

prob = 1-chis_prb(lmlag_rob,1);

result.meth = 'lmlag_rob';
result.lmlag_rob = lmlag_rob;
result.prob = prob;
result.chi1 = 6.635;
result.nobs = n;
result.nvar = k;
