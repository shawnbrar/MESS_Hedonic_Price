% sar_panel_gd3 demo file
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

% calculate mean and std for the effects estimates
% using 2000 MCMC draws

total2 = zeros(ndraw-nomit,n*t);
direct2 = zeros(ndraw-nomit,n*t);
indirect2 = zeros(ndraw-nomit,n*t);


for iter=1:ndraw-nomit
rho = result1.pdraw(iter,1);
beta = result1.bdraw(iter,1);

V = (speye(n*t) - rho*Wtime)\speye(n*t);
bmat  = eye(n*t)*beta(1);
total2(iter,:) = (sum(V*bmat,2))';
direct2(iter,:) = diag(V*bmat);
indirect2(iter,:) = (total2(iter,:) - direct2(iter,:))';
end

direct_int = plims(direct2);
indirect_int = plims(indirect2);
% % plims returns an n*t x 5 matrix with
% % p quantiles from columns of direct2
% % p =[0.005, 0.025, 0.5, 0.975, 0.995];   
% 
tt=1:2*n;
subplot(2,1,1),
plot(tt,direct(tt,1),'.b',tt,direct_int(tt,3),'.r',tt,direct_int(tt,2),'-g',tt,direct_int(tt,4),'-g');
legend('true','estimate','upper0.95','lower0.05');
xlabel('observation-level direct effects');
subplot(2,1,2),
plot(tt,indirect(tt,1),'.b',tt,indirect_int(tt,3),'.r',tt,indirect_int(tt,2),'-g',tt,indirect_int(tt,4),'-g');
legend('true','estimate','upper0.95','lower0.05');
xlabel('observation-level indirect effects');

