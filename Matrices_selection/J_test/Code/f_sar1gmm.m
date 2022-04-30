function f=f_sar1gmm(para_ini,y,x,Wy,Q,P1,We)
n=length(y);
lam=para_ini(1);
bet=para_ini(2:end);
eps=y-lam*Wy-x*bet;
P1eps=P1*eps;
m=[eps'*P1eps; Q'*eps];
f=m'*We*m;
