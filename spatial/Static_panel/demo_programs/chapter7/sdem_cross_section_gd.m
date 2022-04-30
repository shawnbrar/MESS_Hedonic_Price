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
W = make_neighborsw(xc,yc,5); % 5 nearest neighbors W-matrix

Wbig = kron(eye(t),W);

k=4; % 4 explanatory variables
x = [randn(n*t,k)];
beta = [1
        1
        1
        1];
    
theta = 0.5*beta;
    
sige = 1;
rho = 0.8;

u = (speye(n*t) - rho*kron(eye(t),W))\randn(n*t,1)*sqrt(sige);
y = (ones(n,1) + x*beta + Wbig*x*theta + u);

ndraw = 2500;
nomit = 500;
prior.model = 0; % no fixed effects

% Estimation of Bayesian model averaging estimates using three matrices
result = sdem_panel_FE_g(y,[ones(n,1) x],W,t,ndraw, nomit, prior);

vnames = strvcat('y','constant','x1','x2','x3','x4');
prt_panel(result, vnames);

