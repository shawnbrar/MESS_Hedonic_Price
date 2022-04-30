% sar BMA program for sar_conv_panel_bma_gd.m 
clear all;
sd = 221010;
rng(sd);

% estimate all possible models
% with two or more W-matrices
nweights = 5;

% np = 26
n = 400; 
t = 5;
xc = randn(n,1);  % generate 5 W-matrices
yc = randn(n,1);
W1 = make_neighborsw(xc,yc,5); % 5 nearest neighbors W-matrix

xc = randn(n,1);
yc = randn(n,1);
W2 = make_neighborsw(xc,yc,8); % 8 nearest neighbors W-matrix

xc = randn(n,1);
yc = randn(n,1);
W3 = make_neighborsw(xc,yc,10); % 10 nearest neighbors W-matrix

xc = randn(n,1);
yc = randn(n,1);
W4 = make_neighborsw(xc,yc,6); % 6 nearest neighbors W-matrix

xc = randn(n,1);
yc = randn(n,1);
W5 = make_neighborsw(xc,yc,4); % 4 nearest neighbors W-matrix

gamma1 = 0.3; % assign gamma weights
gamma2 = 0.3;
gamma3 = 0.2;
gamma4 = 0.2;
gamma5 = 0.0; % W5 is not really in the model
gtrue = [gamma1
         gamma2
         gamma3
         gamma4
         gamma5];
% 
Wc = gamma1*W1 + gamma2*W2 + gamma3*W3 + gamma4*W4 + gamma5*W5;
u = randn(n,1);

corrcoef([W1*u W2*u W3*u W4*u W5*u])

k=4; % 4 explanatory variables
x = [randn(n*t,k)];
beta = [1
        -1
        -0.5
        1.5];
btrue = beta;        
sige = 1;
strue = sige;
rho = 0.7;
ptrue = rho;
% generate True model
% add fixed effects to the DGP
tts = (1:n)*(1/n);
SFE = kron(ones(t,1),tts');
ttt = (1:t)*(1/t);
TFE = kron(ttt',ones(n,1));

y = (speye(n*t) - rho*kron(eye(t),Wc))\(x*beta + SFE + TFE + randn(n*t,1)*sqrt(sige));

ndraw = 30000;
nomit = 10000;
prior.thin = 5; % retains only 2000 draws from 30,000
                % by skipping every 5
prior.model = 3; % fixed effects

Wmatrices = [W1 W2 W3 W4 W5];

% Estimation of Bayesian model averaging estimates using three matrices
result = sar_conv_panel_bma_g(y,x,Wmatrices,n,t, ndraw, nomit, prior);
vnames = strvcat('y','x1','x2','x3','x4');
prt_panel_bma(result, vnames);

prior.parallel = 0;
result = sar_conv_panel_bma_g(y,x,Wmatrices,n,t, ndraw, nomit, prior);
vnames = strvcat('y','x1','x2','x3','x4');
prt_panel_bma(result, vnames);
