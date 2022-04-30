% sdm BMA program for sdm_conv_panel_bma_gd.m 
clear all;
sd = 221010;
rng(sd);

% estimate all possible models
% with two or more W-matrices
% nweights = 3, so we have 4 models with 2 or more W-matrices

n = 3000; 
t = 1;
m=3;

xc = randn(n,1);  % generate 5 W-matrices
yc = randn(n,1);
W1 = make_neighborsw(xc,yc,5); % 5 nearest neighbors W-matrix

xc = randn(n,1);
yc = randn(n,1);
W2 = make_neighborsw(xc,yc,8); % 8 nearest neighbors W-matrix

xc = randn(n,1);
yc = randn(n,1);
W3 = make_neighborsw(xc,yc,12); % 12 nearest neighbors W-matrix


gamma1 = 0.5; % assign gamma weights
gamma2 = 0.3;
gamma3 = 0.2;
gtrue = [gamma1
         gamma2
         gamma3];
% 
Wc = gamma1*W1 + gamma2*W2 + gamma3*W3;

W1big = kron(eye(t),W1);
W2big = kron(eye(t),W2);
W3big = kron(eye(t),W3);


k=4; % 4 explanatory variables
x = [randn(n*t,k)];
beta = [1
        1
        1
        1];

sige = 1;
rho = 0.4;

u = (speye(n*t) - rho*kron(eye(t),Wc))\randn(n*t,1)*sqrt(sige);
y = (x*beta + u);

ndraw = 50000;
nomit = 10000;
prior.thin = 5; % retains only 8000 draws from 40,000
                % by skipping every 5
prior.model = 0; % no fixed effects

Wmatrices = [W1 W2 W3];

% Estimation of Bayesian model averaging estimates using three matrices
result = sem_conv_panel_g(y,[ones(n,1) x],Wmatrices,n,t, ndraw, nomit, prior);

vnames = strvcat('y','constant','x1','x2','x3','x4');
prt_panel(result, vnames);

