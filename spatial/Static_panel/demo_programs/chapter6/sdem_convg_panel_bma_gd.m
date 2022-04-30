% file: sdem_conv_panel_bma_gd.m
clear all;
sd = 221010;
rng(sd);

% estimate all possible models
% with two or more W-matrices
% nweights = 3, so we have 4 models with 2 or more W-matrices

n = 800; 
t = 5;
xc = randn(n,1);  % generate 5 W-matrices
yc = randn(n,1);
W1 = make_neighborsw(xc,yc,5); % 5 nearest neighbors W-matrix

xc = randn(n,1);
yc = randn(n,1);
W2 = make_neighborsw(xc,yc,8); % 8 nearest neighbors W-matrix

xc = randn(n,1);
yc = randn(n,1);
W3 = make_neighborsw(xc,yc,12); % 12 nearest neighbors W-matrix

m=3;
gamma1 = 0.5; % assign gamma weights
gamma2 = 0.3;
gamma3 = 0.2;
gtrue = [gamma1
         gamma2
         gamma3];
% 
Wc = gamma1*W1 + gamma2*W2 + gamma3*W3;

k=4; % 4 explanatory variables
x = [randn(n*t,k)];
beta = [1
        1
        1
        1];
theta1 = 0.5*beta;
theta2 = -0.75*beta;
theta3 =  beta;

bvec = [beta
        theta1
        theta2
        theta3];

    
% calculate true direct and indirect effects estimates
for ii=1:k
tmp2 = (eye(n)*beta(ii,1) + eye(n)*theta1(ii,1) + eye(n)*theta2(ii,1) + eye(n)*theta3(ii,1));
total_true(ii,1) = mean(sum(tmp2,2));
tmp1 = eye(n)*beta(ii,1); % + eye(n)*theta1(ii,1) + eye(n)*theta2(ii,1) + eye(n)*theta3(ii,1));
direct_true(ii,1) = mean(diag(tmp1));
indirect_true(ii,1) = total_true(ii,1) - direct_true(ii,1);
end

fprintf(1,'true effects estimates \n');
in.cnames = strvcat('direct','indirect','total');
in.rnames = strvcat('variables','x1','x2','x3','x4');

out = [direct_true indirect_true total_true];
mprint(out,in);

sige = 1;
rho = 0.7;
% generate True model
% add fixed effects to the DGP
tts = (1:n)*(1/n);
SFE = kron(ones(t,1),tts');
ttt = (1:t)*(1/t);
TFE = kron(ttt',ones(n,1));

Wmatrices = [kron(eye(t),W1) kron(eye(t),W2) kron(eye(t),W3)];

    Wx = [x kron(speye(t),W1)*x kron(speye(t),W2)*x kron(speye(t),W3)*x];
    Wxb = Wx*bvec;

u = (speye(n*t) - rho*kron(eye(t),Wc))\randn(n*t,1)*sqrt(sige);
y = (Wxb + SFE + TFE + u);

ndraw = 50000;
nomit = 10000;
prior.thin = 5; % retains only 8000 draws from 40,000
                % by skipping every 5
prior.model = 3; % fixed effects


% Estimation of Bayesian model averaging estimates using three matrices
result = sdem_conv_panel_bma_g(y,x,Wmatrices,n,t, ndraw, nomit, prior);

vnames = strvcat('y','x1','x2','x3','x4');
prt_panel_bma(result, vnames);

