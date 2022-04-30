% file: ols_panel_gd.m
clear all;
rng(10203040);

n = 200;
t = 10;

k = 2;
x = randn(n*t,k); % random normal x-variables
beta = ones(k,1); % true beta = 1
sige = 5;         % true noise variance = 5
evec = randn(n*t,1)*sqrt(sige); % random normal disturbances
% fixed effects for regions and time periods
tts = (1:n)*(1/n);
SFE = kron(ones(t,1),tts');
ttt = (1:t)*(1/t);
TFE = kron(ttt',ones(n,1));
% true DGP (data generating process)
y = (x*beta + SFE + TFE + evec);

ndraw = 2500;
nomit = 500;
prior.novi_flag = 1; % homoscedastic model v_{it} = 1
prior.model = 3;     % model with fixed effects for regions and time period
result1 = ols_panel_FE_g(y,x,t,ndraw,nomit,prior);
vnames = strvcat('y','x1','x2');
prt_panel(result1,vnames);

prior2.rval=5;    % heterosedastic model
prior2.model = 1; % model with fixed effects for regions only
% add outliers to y
youtlier = reshape(y,n,t);
youtlier(:,5) = youtlier(:,5) + 10;
youtlier(:,6) = youtlier(:,6) + 10;
yvec = vec(youtlier);
result2 = ols_panel_FE_g(yvec,x,t,ndraw,nomit,prior2);
prt_panel(result2,vnames);
vmean = result2.vmean;
tt=1:n*t;
plot(tt,vmean);
ylabel('v_{it} estimates');
xlabel('n \times t observations');

prior3.model = 1;            % model with fixed effects for regions only
prior3.novi_flag = 1;        % homoscedastic model v_{it} = 1
prior3.beta = ones(2,1)*0.5; % prior means for beta
prior3.bcov = eye(2)*0.001;  % prior variances for beta
prior3.nu = 0.1;             % IG(a,b) a-value uninformative
prior3.d0 = 0.1;             % IG(a,b) b-value uninformative
result3 = ols_panel_FE_g(y,x,t,ndraw,nomit,prior3);
prt_panel(result3,vnames);
