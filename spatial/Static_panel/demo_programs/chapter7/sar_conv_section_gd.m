% sar_conv_panel_gd demo file
clear all;
rng(10203444);

% read an Arcview shape file for 3,111 US counties
map_results = shape_read('../demo_data/uscounties_projected');
latt = map_results.data(:,3);
long = map_results.data(:,4);
n = length(latt);
t = 1;
% plot(long,latt,'.');

[j,W,j] = xy2cont(latt,long); % Delaunay contiguity W  


rho = 0.6;
k = 2;
x = randn(n*t,k);
beta = ones(k,1);
sige = 1;
evec = randn(n*t,1)*sqrt(sige);

gamma1 = 0.3;
gamma2 = 0.7;
gamma3 = 0.0;

y = (speye(n*t) - rho*W)\(ones(n,1)*2 + x*beta + evec);

ndraw = 2000;
nomit = 1000;
prior.model = 0;
prior.thin = 1;
% prior.plt_flag = 1;


result1 = sar_panel_FE_g(y,[ones(n,1) x],W,t,ndraw,nomit,prior);
vnames = strvcat('y','constant','x1','x2');
prt_panel(result1,vnames);



