% sar_conv_panel_gd2 demo file
clear all;
rng(10203040);

map_results = shape_read('../demo_data/usstates49');
latt = [map_results.data(1:8,2)
       map_results.data(10:end,2)]; % skip Washington DC
long = [map_results.data(1:8,3)
       map_results.data(10:end,3)];
       
n = length(latt);
t = 20;

W6 = make_neighborsw(latt,long,6);% 6 nearest neighbors W

Wdist = distance(latt,long) + eye(n); 
Winv_distance = zeros(n,n);
dmax = max(max(Wdist));
for i=1:n
    for j=1:n
        if Wdist(i,j) ~= 0
            Winv_distance(i,j) = 1/Wdist(i,j);
        else
            Winv_distance(i,j) = 1/dmax;
        end
    end
end
% create inverse distance with a 6 neighbor cut-off W
Wtmp = Winv_distance.*W6;
Winv_distance = normw(Wtmp);


% state-to-state commodity flows, 2017
[a,b] = xlsread('../demo_data/cflows_2017.xlsx');
% set main diagonal (intrastate flows) to zero
diaga = diag(a);
W = a - diag(diaga);

Wcom_flows = normw(W); % row-normalize
% eliminate small elements
for i=1:n
    for j=1:n
        if Wcom_flows(i,j) < 0.005
            Wcom_flows(i,j) = 0;
        end
    end
end

Wcom_flows = normw(Wcom_flows);

% miles of common borders between states
[a,b] = xlsread('../demo_data/states_borders.xlsx');
snames = strvcat(b(2:end,:));
% only upper triangle
Wmiles = a(:,2:end);
% make it symmetric
for i=1:48
    for j=1:48
        if Wmiles(i,j) > 0
            Wmiles(j,i) = Wmiles(i,j);
        end
    end
end
Wborder_miles = normw(Wmiles); % row-normalize

% 48 x 48 binary contiguity matrix for states
[a,b] = xlsread('../demo_data/Wcont48.xlsx');
Wcontiguity = normw(a);

% find correlation between W-matrices
u = randn(n,1);
corr = corrcoef([Wcom_flows*u Wborder_miles*u Wcontiguity*u Winv_distance*u]);
inc.cnames = strvcat('Wcom','Wborder','Wcont','Wdist');
inc.rnames = strvcat('Correlation','Wcom','Wborder','Wcont','Wdist');
mprint(corr,inc);

rho = 0.7;
k = 2;
x = randn(n*t,k);
beta = ones(k,1);
sige = 0.1;
evec = randn(n*t,1)*sqrt(sige);

gamma1 = 0.2;
gamma2 = 0.5;
gamma3 = 0.1;
gamma4 = 0.2;

Wc = gamma1*kron(eye(t),Wcom_flows) + gamma2*kron(eye(t),Wborder_miles) + ...
    gamma3*kron(eye(t),Wcontiguity)+ gamma4*kron(eye(t),Winv_distance);

% add fixed effects to the DGP
tts = (1:n)*(1/n);
SFE = kron(ones(t,1),tts');
ttt = (1:t)*(1/t);
TFE = kron(ttt',ones(n,1));

y = (speye(n*t) - rho*Wc)\(x*beta + SFE + TFE + evec);

ndraw = 20000;
nomit = 2000;
prior.model = 3;
prior.novi_flag = 1;
prior.thin = 4;
prior.plt_flag = 1;

Wtrue = gamma1*Wcom_flows + gamma2*Wborder_miles + ...
    gamma3*Wcontiguity+ gamma4*Winv_distance;

result1 = sar_panel_FE_g(y,x,Wtrue,t,ndraw,nomit,prior);
vnames = strvcat('y','x1','x2');
fprintf(1,'Estimates based on true Wc \n');
prt_panel(result1,vnames);

Wmatrices = [kron(eye(t),Wcom_flows) kron(eye(t),Wborder_miles) kron(eye(t),Wcontiguity) kron(eye(t),Winv_distance)];

result2 = sar_conv_panel_g(y,x,Wmatrices,n,t,ndraw,nomit,prior);
fprintf(1,'Estimates based on 4 W-matrices \n');
prt_panel(result2,vnames);

Wmatrices = [kron(eye(t),Wcom_flows) kron(eye(t),Wborder_miles)];

result3 = sar_conv_panel_g(y,x,Wmatrices,n,t,ndraw,nomit,prior);
fprintf(1,'Estimates based on 2 W-matrices \n');
prt_panel(result3,vnames);

