% sar_conv_cross_section_gd demo file
clear all;
rng(10203444);

% read an Arcview shape file for 3,111 US counties
map_results = shape_read('../demo_data/uscounties_projected');
latt = map_results.data(:,3);
long = map_results.data(:,4);
n = length(latt);
t = 1;
% plot(long,latt,'.');

[j,Wcont,j] = xy2cont(latt,long); % Delaunay contiguity W  

W6 = make_neighborsw(latt,long,6);% 6 nearest neighbors W

Wdist = distance(latt,long) + eye(n);
% create inverse distance with a 6 neighbor cut-off W
Wcut = (ones(n,n)./Wdist).*W6;
Wdist = normw(Wcut);

rho = 0.6;
k = 2;
x = randn(n*t,k);
beta = ones(k,1);
sige = 1;
evec = randn(n*t,1)*sqrt(sige);

gamma1 = 0.3;
gamma2 = 0.7;
gamma3 = 0.0;

Wc = gamma1*kron(eye(t),Wcont) + gamma2*kron(eye(t),W6) + gamma3*kron(eye(t),Wdist);

y = (speye(n*t) - rho*Wc)\(x*beta + evec);

ndraw = 20000;
nomit = 10000;
prior.model = 0;
prior.thin = 5;
% prior.plt_flag = 1;


Wmatrices = [Wcont W6 Wdist];

result1 = sar_conv_panel_g(y,x,Wmatrices,n,t,ndraw,nomit,prior);
vnames = strvcat('y','x1','x2');
prt_panel(result1,vnames);


result2 = sar_conv_panel_bma_g(y,x,Wmatrices,n,t,ndraw,nomit,prior);
vnames = strvcat('y','x1','x2');
prt_panel_bma(result2,vnames);

