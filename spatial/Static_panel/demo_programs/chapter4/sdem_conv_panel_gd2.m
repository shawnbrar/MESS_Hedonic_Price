% sdem_panel_gd2 demo file
clear all;
rng(30203040);

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
Wc = gamma1*kron(eye(t),W1) + gamma2*kron(eye(t),W2) + gamma3*kron(eye(t),W3);

k=4; % 4 explanatory variables
x = [randn(n*t,k)];
beta = ones(k,1);
theta1 = -0.75*ones(k,1);
theta2 = -0.25*ones(k,1);
theta3 = -1*ones(k,1);

bvec = [beta
        theta1
        theta2
        theta3];

rho = 0.7;
% generate True model
% add fixed effects to the DGP
tts = (1:n)*(1/n);
SFE = kron(ones(t,1),tts');
ttt = (1:t)*(1/t);
TFE = kron(ttt',ones(n,1));

Wmatrices = [kron(eye(t),W1) kron(eye(t),W2) kron(eye(t),W3)];

% the DGP must match the order in which
% the function sdm_conv_panel_g()
% processes the x, W*x variables

% the DGP x-variables must match the order in which
% the function sdm_conv_panel_g()
% processes the x, W*x variables
Wx = x;
begi = 1;
endi = n*t;
for ii=1:m
        Wx = [Wx Wmatrices(:,begi:endi)*x];
    begi = begi + n*t;
    endi = endi + n*t;
end

Wxb = Wx*bvec;


Wc = gamma1*kron(eye(t),W1) + gamma2*kron(eye(t),W2) + gamma3*kron(eye(t),W3);

sige = 1;

evec = (speye(n*t) - rho*Wc)\(sqrt(sige)*randn(n*t,1));

y = (Wxb + SFE + TFE + evec);


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

prior.model = 3;
prior.plt_flag = 0;

ndraw = 20000;
nomit = 10000;
prior.thin = 2; 
% we use only 4,000 of the 20,000 retained draws
% to calculate effects estimates and posterior means, std

Wmatrices = [W1 W2 W3];

result1 = sdem_conv_panel_g(y,x,Wmatrices,n,t,ndraw,nomit,prior);
vnames = strvcat('y','x1','x2','x3','x4');
prt_panel(result1,vnames);


