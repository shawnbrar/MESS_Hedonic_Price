% sdm_conv_cross_section_gd demo file
% demonstrate using the function to produce
% estimates for a cross-sectional model

clear all;
rng(30203040);

n = 3000;
t = 1;

rho = 0.2;
k = 2;
x = [randn(n*t,k)];
beta = ones(k,1);
theta1 = -0.75*ones(k,1);
theta2 = -0.25*ones(k,1);

bvec = [beta
        theta1
        theta2];

sige = 1;
evec = randn(n*t,1)*sqrt(sige);
latt = rand(n,1);
long = rand(n,1);
W1 = make_neighborsw(latt,long,2);
latt = rand(n,1); % A different set of latt-long
long = rand(n,1); % coordinates
W2 = make_neighborsw(latt,long,6);
latt = rand(n,1); % A different set of latt-long
long = rand(n,1); % coordinates
W3 = make_neighborsw(latt,long,12);

gamma1 = 0.2;
gamma2 = 0.8;

    Wx = [x kron(speye(t),W1)*x kron(speye(t),W2)*x];
    Wxb = Wx*bvec;

% calculate true direct and indirect effects estimates
Wc = gamma1*W1 + gamma2*W2;

direct_true = zeros(k,1);
indirect_true = zeros(k,1);
total_true = zeros(k,1);

B = (speye(n) - rho*Wc)\(speye(n));

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
            
y = (speye(n*t) - rho*Wc)\(ones(n,1)*2 + Wxb + evec);

prior.model = 0;
prior.plt_flag = 0;
ndraw = 20000;
nomit = 10000;
prior.thin = 2;

Wmatrices = [W1 W2 W3];

result1 = sdm_conv_panel_g(y,[ones(n,1) x],Wmatrices,n,t,ndraw,nomit,prior);
vnames = strvcat('y','constant','x1','x2');
prt_panel(result1,vnames);

result2 = sdm_conv_panel_bma_g(y,[ones(n,1) x],Wmatrices,n,t,ndraw,nomit,prior);
vnames = strvcat('y','constant','x1','x2');
prt_panel_bma(result2,vnames);
