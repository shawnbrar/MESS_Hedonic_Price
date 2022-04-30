% sar_conv_panel_gd demo file
clear all;
rng(10203444);

% read an Arcview shape file for 3,111 US counties
map_results = shape_read('../demo_data/uscounties_projected');
latt = map_results.data(:,3);
long = map_results.data(:,4);
n = length(latt);
t = 20;

[j,Wcont,j] = xy2cont(latt,long); % Delaunay contiguity W  

W6 = make_neighborsw(latt,long,6);% 6 nearest neighbors W

Wdist = distance(latt,long) + eye(n);
% create inverse distance with a 6 neighbor cut-off W
Wcut = (ones(n,n)./Wdist).*W6;
Wdist = normw(Wcut);

rho = 0.1;
k = 2;
x = randn(n*t,k);
beta = ones(k,1);
sige = 1;
evec = randn(n*t,1)*sqrt(sige);

gamma1 = 0.3;
gamma2 = 0.6;
gamma3 = 0.1;

Wc = gamma1*kron(eye(t),Wcont) + gamma2*kron(eye(t),W6) + gamma3*kron(eye(t),Wdist);

% add fixed effects to the DGP
tts = (1:n)*(1/n);
SFE = kron(ones(t,1),tts');
ttt = (1:t)*(1/t);
TFE = kron(ttt',ones(n,1));

y = (speye(n*t) - rho*Wc)\(x*beta + SFE + TFE + evec);

ndraw = 50000;
nomit = 20000;
prior.model = 3;
prior.novi_flag = 1;
prior.thin = 10;
% prior.plt_flag = 1;

Wmatrices = [kron(eye(t),Wcont) kron(eye(t),W6) kron(eye(t),Wdist)];

result1 = sar_conv_panel_g(y,x,Wmatrices,n,t,ndraw,nomit,prior);
vnames = strvcat('y','x1','x2');
fprintf(1,'sar_conv_panel_g with rho = 0.1 \n');
prt_panel(result1,vnames);


rho = -0.2;
k = 2;
x = randn(n*t,k);
beta = ones(k,1);
sige = 1;
evec = randn(n*t,1)*sqrt(sige);

gamma1 = 0.3;
gamma2 = 0.6;
gamma3 = 0.1;

Wc = gamma1*kron(eye(t),Wcont) + gamma2*kron(eye(t),W6) + gamma3*kron(eye(t),Wdist);

% add fixed effects to the DGP
tts = (1:n)*(1/n);
SFE = kron(ones(t,1),tts');
ttt = (1:t)*(1/t);
TFE = kron(ttt',ones(n,1));

y = (speye(n*t) - rho*Wc)\(x*beta + SFE + TFE + evec);

Wmatrices = [kron(eye(t),Wcont) kron(eye(t),W6) kron(eye(t),Wdist)];

result2 = sar_conv_panel_g(y,x,Wmatrices,n,t,ndraw,nomit,prior);
vnames = strvcat('y','x1','x2');
fprintf(1,'sar_conv_panel_g with rho = -0.2 \n');
prt_panel(result2,vnames);


