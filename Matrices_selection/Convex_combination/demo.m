%%%%%%%%%%%%
% Demo file%
%%%%%%%%%%%%

clear all;
%% DGP
% set seed
S=123;
rng(S);
n = 2000;
In=eye(n);
x=[ones(n,1), randn(n,2)];
bet=[1;2;3];
xb=x*bet;
rh=0.5;
g1=0.2;
g2=0.5;
g3=1-g1-g2;

W1=make_neighborsw(rand(n,1),rand(n,1),5);
W2=make_neighborsw(rand(n,1),rand(n,1),5);
W3=make_neighborsw(rand(n,1),rand(n,1),5);
Wc=g1*W1 +g2*W2 + g3*W3;

A= In-rh*Wc;
Ai=In/A;
e=randn(n,1)*0.4;
y = Ai*(xb + e);

%% Model estimation 
vnames=strvcat('y','constant','x_1','x_2');
ndraw=25000;
nomit=5000;
prior.plt=1;% Plot MCMC draws

% Estimation of a single convex combination made of the three matrices 
Wmatrices=[W1 W2 W3];
res_3=sar_conv_g(y,x,Wmatrices,ndraw, nomit,prior);
prt_sar_conv_g(res_3,vnames);

% Estimation of a single convex combination made of the first 2 matrices 
Wmatrices=[W1 W2];
res_2=sar_conv_g(y,x,Wmatrices,ndraw, nomit,prior);
prt_sar_conv_g(res_2,vnames);


% Estimation of a conventional SAR model with W1

res_1=sar_g(y,x,W1,ndraw, nomit,prior);
prt_sar_g(res_1,vnames);

% Estimation of Bayesian model averaging estimates using three matrices
res_bma=sar_conv_bma_g(y,x,Wmatrices, ndraw, nomit, prior);
prt_sar_conv_bma_g(res_bma, vnames);