% sar_panel_g demo file

clear all;

rng(10203040);

n = 200;
t = 10;

rho = 0.6;
k = 2;
x = randn(n*t,k);
beta = ones(k,1);
sige = 1;
evec = randn(n*t,1)*sqrt(sige);

latt = rand(n,1);
long = rand(n,1);

W = make_neighborsw(latt,long,5);

Wbig = kron(eye(t),W);

% add fixed effects to the DGP

tts = (1:n)*(1/n);
SFE = kron(ones(t,1),tts');
ttt = (1:t)*(1/t);
TFE = kron(ttt',ones(n,1));

y = (speye(n*t) - rho*Wbig)\(x*beta + SFE + TFE + evec);

prior.model = 3;
result1 = sar_panel_FE(y,x,W,t,prior);
vnames = strvcat('y','x1','x2');
prt_panel(result1,vnames);

ndraw = 2500;
nomit = 500;
prior2.novi_flag = 1;
prior2.model = 3;
result2 = sar_panel_FE_g(y,x,W,t,ndraw,nomit,prior2);
prt_panel(result2,vnames);

result3 = sar_panel_FE_g(y,x,Wbig,t,ndraw,nomit,prior2);
prt_panel(result3,vnames);


