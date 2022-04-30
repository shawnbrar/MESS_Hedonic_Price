% sdem_panel_g demo file

clear all;
rng(10203040);

n = 100;
t = 20;
rho = 0.7;
k = 2;
x = randn(n*t,k);
beta = ones(k,1);
beta = [beta
        -ones(k,1)];
sige = 0.1;
evec = randn(n*t,1)*sqrt(sige);

latt = rand(n,1);
long = rand(n,1);

W = make_neighborsw(latt,long,5);

Wbig = kron(eye(t),W);

xo = x;
x = [x Wbig*x];

% add fixed effects to the DGP

tts = (1:n)*(1/n);
SFE = kron(ones(t,1),tts');
ttt = (1:t)*(1/t);
TFE = kron(ttt',ones(n,1));

u = (speye(n*t) - rho*Wbig)\evec;
y = (x*beta + SFE + TFE + u);

info.model = 3;
result1 = sdem_panel_FE(y,xo,W,t,info);
vnames = strvcat('y','x1','x2');
prt_panel(result1,vnames);

ndraw = 2500;
nomit = 500;
prior.novi_flag = 1;
prior.model = 3;
result2 = sdem_panel_FE_g(y,xo,W,t,ndraw,nomit,prior);
prt_panel(result2,vnames);

prior2.rval = 5;
prior2.model = 3;
result3 = sdem_panel_FE_g(y,xo,W,t,ndraw,nomit,prior2);
prt_panel(result3,vnames);

prior3.novi_flag = 1;
prior3.model = 3;
k = size(x,2);
prior3.beta = zeros(k,1);
prior3.bcov = eye(k)*10;
result4 = sdem_panel_FE_g(y,xo,W,t,ndraw,nomit,prior3);
prt_panel(result4,vnames);

prior4.model = 3;
prior4.beta = ones(4,1);
prior4.bcov = eye(4)*0.001;
result5 = sdem_panel_FE_g(y,xo,W,t,ndraw,nomit,prior4);
prt_panel(result5,vnames);



