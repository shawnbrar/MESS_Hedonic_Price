% sdem BMA program for sdem_conv_cross_section_gd.m 
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

k=2; % 4 explanatory variables
x = [randn(n*t,k)];
beta = [1
        1];
theta1 = 0.5*beta;
theta2 = 1*beta;
theta3 = -1*beta;

bvec = [beta
        theta1
        theta2
        theta3];


sige = 1;
rho = 0.6;

    
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
in.rnames = strvcat('variables','x1','x2');

out = [direct_true indirect_true total_true];
mprint(out,in);

    Wx = [x kron(speye(t),W1)*x kron(speye(t),W2)*x kron(speye(t),W3)*x];
    Wxb = Wx*bvec;


u = (speye(n*t) - rho*kron(eye(t),Wc))\randn(n*t,1)*sqrt(sige);
y = (ones(n,1)*2.0 + Wxb + u);

ndraw = 20000;
nomit = 10000;
prior.thin = 5; % retains only 2000 draws from 10,000
                % by skipping every 5
prior.model = 0; % no fixed effects
% prior.plt_flag = 1;

Wmatrices = [W1 W2 W3];

% Estimation of Bayesian model averaging estimates using three matrices
result = sdem_conv_panel_g(y,[ones(n,1) x],Wmatrices,n,t, ndraw, nomit, prior);
vnames = strvcat('y','constant','x1','x2');
prt_panel(result, vnames);

result = sdem_conv_panel_bma_g(y,[ones(n,1) x],Wmatrices,n,t, ndraw, nomit, prior);
vnames = strvcat('y','constant','x1','x2');
prt_panel_bma(result, vnames);
