% sdm_panel_gd demo file
% demonstrate using the function to produce
% estimates for a cross-sectional model

clear all;
rng(30203040);

n = 1000;
t = 1;

rho = 0.5;
k = 2;
x = [randn(n*t,k)];
beta = ones(k,1);
theta1 = -0.5*ones(k,1);

sige = 0.5;
evec = randn(n*t,1)*sqrt(sige);
latt = rand(n,1);
long = rand(n,1);

W = make_neighborsw(latt,long,6);


Wbig = kron(eye(t),W);

y = (speye(n*t) - rho*W)\(ones(n,1)*2 + x*beta + Wbig*x*theta1 + evec);

% calculate true direct and indirect effects estimates
direct_true = zeros(k,1);
indirect_true = zeros(k,1);
total_true = zeros(k,1);

B = (speye(n) - rho*W);

for ii=1:k
tmp2 = B\(eye(n)*beta(ii,1) + W*theta1(ii,1));
total_true(ii,1) = mean(sum(tmp2,2));
tmp1 = B\(eye(n)*beta(ii,1) + W*theta1(ii,1));
direct_true(ii,1) = mean(diag(tmp1));
indirect_true(ii,1) = total_true(ii,1) - direct_true(ii,1);
end

fprintf(1,'true effects estimates \n');
in.cnames = strvcat('direct','indirect','total');
in.rnames = strvcat('variables','x1','x2');

out = [direct_true indirect_true total_true];
mprint(out,in);


prior.model = 0;

ndraw = 20000;
nomit = 10000;


result1 = sdm_panel_FE_g(y,[ones(n,1) x],W,t,ndraw,nomit,prior);
vnames = strvcat('y','constant','x1','x2');
prt_panel(result1,vnames);

