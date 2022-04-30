% Example with FD data
% This dataset includes data on Investment and saving for 24 countries over
% the period 1960 to 2000. 
%  We wish to test which spatial model should be considered. 
clear 
load Example_FH.mat

%  Estimation on subsample
id=d6070t(:,26);
vnamesw60=strvcat('investmentw60', 'savingsw60');
resw60 = pfixed(d6070t(:,1:2),id);
prt_c_panel(resw60,vnamesw60);

% % From 1971 to 1985 (15 periods)
% 
id=d7185t(:,26);
vnamesw71=strvcat('investmentw71', 'savingsw71');
resw71 = pfixed(d7185t(:,1:2),id);
prt_c_panel(resw71,vnamesw71);
% 
% % From 1986 to 2000 (15 periods)
% 
id=d8600t(:,26);
vnamesw86=strvcat('investmentw86', 'savingsw86');
resw86 = pfixed(d8600t(:,1:2),id);
prt_c_panel(resw86,vnamesw86);


% Definition of spatial weight matrices
N=length(lat);
Wdia1=dinvarc(lat,long,1);
Wdia1=normw(Wdia1);
% Wdia2=dinvarc(lat,long,2);
Wdeia=dexparc(lat,long,1);
Wdeia=normw(Wdeia);
% l=eig(Wdia1);
% t=max(abs(l));
% Wdt=Wdia1/t;
%Creation of the k-nearest neighbours.
m=7; %Number of nearest neighbours


W = make_neighborsw(lat,long,m);
W=normw(W);
y = d8600t(:,1); % Investment
x = d8600t(:,2); % Savings

% Estimation of spatial models
info.model=1; % individual fixed effects.
info.lflag=0; 
T=length(y)/N;
results=sar_panel_FE(y, x,W,T,info); 
vnames=strvcat('investment','savings');
prt_panel(results,vnames,1);

Wsavings=kron(eye(T),W)*x;
Z=[x Wsavings];
%% Spatial autocorrelation tests
T=length(y)/size(W,1);
res=sar_panel_FE(y,x,W,T,info);
vnamesw86=strvcat('investmentw86', 'savingsw86');
prt_panel(res,vnamesw86)
res1=sem_panel_FE(y,x,W,T,info);
prt_panel(res1,vnamesw86)


results1   =  lm_f_joint(y,x,W,W,N);
prt(results1)
results2  =  lm_f_err(y,x,W,N);
prt(results2)
results3  =  lm_f_sar(y,x,W,N);
prt(results3)
results4  =lm_f_sar_c(y,x,W,W,N);
prt(results4)
results5  =lm_f_err_c(y,x,W,W,N);
prt(results5)









