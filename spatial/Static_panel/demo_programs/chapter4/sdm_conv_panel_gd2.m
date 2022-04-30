% sdm_conv_panel_gd2 demo file
clear all;
rng(10203040);

n = 500;
t = 5;

rho = 0.7;
k = 2;
x = randn(n*t,k);
beta = [1
        0.5];
theta1 = -0.75*ones(k,1);
theta2 = 0.5*ones(k,1);

% (1/(1-rho))*(beta + theta1 + theta2)

bvec = [beta
        theta1
        theta2];

latt = rand(n,1);
long = rand(n,1);

W1 = make_neighborsw(latt,long,2);

latt = rand(n,1); % A different set of latt-long
long = rand(n,1); % coordinates

W2 = make_neighborsw(latt,long,6);

m = 2;

gamma1 = 0.2;
gamma2 = 0.8;

gamma = [gamma1
         gamma2];

Wc = gamma1*kron(eye(t),W1) + gamma2*kron(eye(t),W2);

% add fixed effects to the DGP
tts = (1:n)*(1/n);
SFE = kron(ones(t,1),tts');
ttt = (1:t)*(1/t);
TFE = kron(ttt',ones(n,1));

Wmatrices = [kron(eye(t),W1) kron(eye(t),W2)];

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

sige = 1;
evec = sqrt(sige)*randn(n*t,1);

y = (speye(n*t) - rho*Wc)\(Wxb + SFE + TFE + evec);

% calculate true direct and indirect effects estimates
Wc_small = gamma1*W1 + gamma2*W2;

direct_true = zeros(k,1);
indirect_true = zeros(k,1);
total_true = zeros(k,1);

B = (speye(n) - rho*Wc_small)\(speye(n));

for ii=1:k
tmp2 = B*(eye(n)*beta(ii,1) + W1*theta1(ii,1) + W2*theta2(ii,1));
total_true(ii,1) = mean(sum(tmp2,2));
tmp1 = B*(eye(n)*beta(ii,1) + W1*theta1(ii,1) + W2*theta2(ii,1));
direct_true(ii,1) = mean(diag(tmp1));
indirect_true(ii,1) = total_true(ii,1) - direct_true(ii,1);
end

fprintf(1,'true effects estimates \n');
in.cnames = strvcat('direct','indirect','total');
in.rnames = strvcat('variables','x1','x2');

out = [direct_true indirect_true total_true];
mprint(out,in);

Wmatrices = [W1 W2];
ndraw = 20000;
nomit = 10000;
prior.model=3;
prior.plt_flag = 0;
prior.thin = 5;
result1 = sdm_conv_panel_g(y,x,Wmatrices,n,t,ndraw,nomit,prior);
vnames = strvcat('y','x1','x2');
prt_panel(result1,vnames);

bdraw = result1.bdraw;
pdraw = result1.pdraw;
gdraw = result1.gdraw;

bmean = mean(bdraw);
[bmean' bvec]


