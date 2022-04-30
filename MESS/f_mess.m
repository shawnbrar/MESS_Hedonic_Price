function llik=f_mess(parm,y,x,W)
% Purpose = concentrated function to minimize over in order to get an
% estimated value of alpha

% y is the vector of the dependent variable (N x 1)
% x is the  matrix of independent variables (N x k)
% W is the weight matrix

[n,k]=size(x);
In=eye(n);
xpx=x'*x;
xpxxp=inv(xpx)*x';
H = In-x*xpxxp;
e=expm(parm*W)*y;
ephe=e'*H*e;
llik = n/2*log(ephe);
