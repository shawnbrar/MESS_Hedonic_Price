% sar_panel_gd2 demo file
clear all;
rng(10203444);

n = 200;
t = 6;

rho = 0.7;
k = 2;
x = randn(n*t,k);
beta = ones(k,1);
sige = 1;
evec = randn(n*t,1)*sqrt(sige);

latt = rand(n,1);
long = rand(n,1);

W1 = make_neighborsw(latt,long,10);
W2 = make_neighborsw(latt,long,3);

Wtime = blockdiag(W1,W2,W1,W2,W1,W2);

% add fixed effects to the DGP
tts = (1:n)*(1/n);
SFE = kron(ones(t,1),tts');
ttt = (1:t)*(1/t);
TFE = kron(ttt',ones(n,1));

y = (speye(n*t) - rho*Wtime)\(x*beta + SFE + TFE + evec);

% calculate true observation-level effects estimates
V = (speye(n*t) - rho*Wtime)\speye(n*t);
bmat  = eye(n*t)*beta(1);
total = sum(V*bmat,2);
direct = diag(V*bmat);
indirect = total - direct;
direct_mean = mean(direct);
indirect_mean = mean(indirect);
total_mean = mean(total);

[direct_mean indirect_mean total_mean]

ndraw = 2500;
nomit = 500;
info.novi_flag = 1;
info.model = 3;

result1 = sar_panel_FE_g(y,x,Wtime,t,ndraw,nomit,info);
prt_panel(result1);

rho = result1.rho;
beta = result1.beta(1,1);

V = (speye(n*t) - rho*Wtime)\speye(n*t);
bmat  = eye(n*t)*beta(1);
total2 = sum(V*bmat,2);
direct2 = diag(V*bmat);
indirect2 = total2 - direct2;
direct_mean = mean(direct2);
indirect_mean = mean(indirect2);
total_mean = mean(total2);

tt=1:n*t;
subplot(2,1,1),
plot(tt,direct,'.b',tt,direct2,'.r');
legend('true','estimate');
xlabel('observation-level direct effects');
subplot(2,1,2),
plot(tt,indirect,'.b',tt,indirect2,'.r');
legend('true','estimate');
xlabel('observation-level indirect effects');
