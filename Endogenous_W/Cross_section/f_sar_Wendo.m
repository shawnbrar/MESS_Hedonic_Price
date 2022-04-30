function llike = f_sar_Wendo(theta,detval,n,Zn,x2,y,W,x1)
% PURPOSE: evaluates concentrated log-likelihood for the 
%  spatial autoregressive model using sparse matrix algorithms
% ---------------------------------------------------
%  USAGE:llike = f_sar(rho,detval,epe0,eped,epe0d,n)
%  where: rho  = spatial autoregressive parameter
%         detval = a matrix with vectorized log-determinant information
%         epe0   = see below
%         eped   = see below
%         eoe0d  = see below
%         n      = # of obs
%          b0 = AI*xs'*ys;
%          bd = AI*xs'*Wys;
%          e0 = ys - xs*b0;
%          ed = Wys - xs*bd;
%          epe0 = e0'*e0;
%          eped = ed'*ed;
%          epe0d = ed'*e0;
% ---------------------------------------------------
%  RETURNS: a  scalar equal to minus the log-likelihood
%           function value at the parameter rho
% ---------------------------------------------------                         
%  NOTE: this is really two functions depending
%        on nargin = 3 or nargin = 4 (see the function)
%  --------------------------------------------------
%  SEE ALSO: sar, f_far, f_sac, f_sem
% ---------------------------------------------------

% written by: James P. LeSage 1/2000
% University of Toledo
% Department of Economics
% Toledo, OH 43606
% jlesage@spatial-econometrics.com



rho=theta(1,1);
beta=theta(2:size(x1,2)+1,1);
Gamma=theta(size(x1,2)+2:size(x1,2)+1+size(x2,2),1); %hypothèse que p2=1
sig2_ep=theta(size(x1,2)+2+size(x2,2),1); %hypothèse que p2=1
sig2_v=theta(size(x1,2)+2+size(x2,2)+1,1);
coeffcorr=theta(size(x1,2)+2+size(x2,2)+2,1); %hypothèse que p2=1



if nargin == 8
% gsize = detval(2,1) - detval(1,1);
% % Note these are actually log detvalues
% i1 = find(detval(:,1) <= rho + gsize);
% i2 = find(detval(:,1) <= rho - gsize);
% i1 = max(i1);
% i2 = max(i2);
% index = round((i1+i2)/2);
% if isempty(index)
% index = 1;
% end;
% detm = detval(index,2);
%
w=eig(W);
sum1=0;
for i =1:n
    p=1-rho*w(i);
    sum1=sum1+log(p);
end

ep=(Zn-x2*Gamma);
Sn=(eye(n) - rho*W);
eta=Sn*y-x1*beta-ep*coeffcorr*sqrt(sig2_v)/sqrt(sig2_ep);

%fonction de vraissemblance concentré:
%on a ici -ln(L) car on utilise une fonction de minimisation, or il faut maximiser la log-vraissemblance ln(L)
%llike = n/2*log(sig2_ep*sig2_v*(1-coeffcorr^2))  - detm + 1/(2*sig2_ep)*(ep')*ep + 1/(2*sig2_v*(1-coeffcorr^2))*(eta')*eta;
llike = n/2*log(sig2_ep*sig2_v*(1-coeffcorr^2))  - sum1 + 1/(2*sig2_ep)*(ep')*ep + 1/(2*sig2_v*(1-coeffcorr^2))*(eta')*eta;





else

error('f_sar_Wendo: Wrong # of input arguments');

end;
