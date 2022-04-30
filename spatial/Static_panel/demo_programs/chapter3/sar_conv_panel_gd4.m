% sar_conv_panel_gd demo file
clear all;
rng(10203444);

% read an Arcview shape file for 3,111 US counties
map_results = shape_read('../demo_data/uscounties_projected');
latt = map_results.data(:,3);
long = map_results.data(:,4);
n = length(latt);
t = 20;

W6 = make_neighborsw(latt,long,6);% 6 nearest neighbors W
W7 = make_neighborsw(latt,long,7);% 7 nearest neighbors W
W20 = make_neighborsw(latt,long,20);% 20 nearest neighbors W

u = randn(n,1);
corr= corrcoef([W6*u W7*u W20*u]);
    
% find correlation between W-matrices
inc.cnames = strvcat('W6','W7','W20');
inc.rnames = strvcat('Correlation','W6','W7','W20');
mprint(corr,inc);


rho = 0.6;
k = 2;
x = randn(n*t,k);
beta = ones(k,1);
sige = 100;
evec = randn(n*t,1)*sqrt(sige);

gamma1 = 0.3;
gamma2 = 0.7;

Wc = gamma1*kron(eye(t),W6) + gamma2*kron(eye(t),W7);

% add fixed effects to the DGP
tts = (1:n)*(1/n);
SFE = kron(ones(t,1),tts');
ttt = (1:t)*(1/t);
TFE = kron(ttt',ones(n,1));

y = (speye(n*t) - rho*Wc)\(x*beta + SFE + TFE + evec);

ndraw = 20000;
nomit = 10000;
prior.model = 3;
prior.novi_flag = 1;
prior.thin = 4;
% prior.plt_flag = 1;

Wmatrices = [kron(eye(t),W6) kron(eye(t),W7)];

result1 = sar_conv_panel_g(y,x,Wmatrices,n,t,ndraw,nomit,prior);
vnames = strvcat('y','x1','x2');
fprintf(1,'sar_conv_panel_g with W6, W7 \n');
prt_panel(result1,vnames);

Wc = gamma1*kron(eye(t),W6) + gamma2*kron(eye(t),W20);
y = (speye(n*t) - rho*Wc)\(x*beta + SFE + TFE + evec);

Wmatrices = [kron(eye(t),W6) kron(eye(t),W20)];

result2 = sar_conv_panel_g(y,x,Wmatrices,n,t,ndraw,nomit,prior);
vnames = strvcat('y','x1','x2');
fprintf(1,'sar_conv_panel_g with W6, W20 \n');
prt_panel(result2,vnames);


