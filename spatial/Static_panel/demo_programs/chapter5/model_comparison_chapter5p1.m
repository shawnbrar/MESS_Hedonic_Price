% model_comparison_chapter5p1.m demo file
clear all;
[cclaims,b] = xlsread('../demo_data/weekly.xlsx',1);
% read data from sheet 1 of Excel spreadsheet
% growth rate of unemployment 2019-2020 from same week, previous year
snames = strvcat(b(2:end,1)); % 48 state names
tnames = strvcat(b(1,2:end)); % 51 week labels
[N,T] = size(cclaims);
[jobposts,b] = xlsread('../demo_data/weekly.xlsx',2);
% read data from sheet 2 of Excel spreadsheet
% change in job offers from 1st week of 2020
[athome,b] = xlsread('../demo_data/weekly.xlsx',3);
% read data from sheet 3 of Excel spreadsheet
% growth rate of percent population at home
% 2019-2020 from same week, previous year
[a,b] = xlsread('../demo_data/Wcont48.xlsx');
% 48 x 48 contiguity matrix for states
W = normw(a);

y = vec(cclaims);
x = [vec(jobposts) vec(athome)];


model = 0; % no fixed effects
[ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,x,N,T,model);

info.lflag = 0; % exact log-determinant
% info.lflag = 1; uses Pace and Barry approximation
% which is faster for large data samples

result1 = lmarginal_static_panel(ywith,xwith,W,N,T,info); 

fprintf(1,'no fixed effects: log marginal likelihoods and model probabilities \n');
in.cnames = strvcat('log-marginal','model probs');
in.rnames = strvcat('model','slx','sdm','sdem');
in.width = 10000;
in.fmt = '%10.4f';
out1 = [result1.lmarginal result1.probs];
mprint(out1,in);

model = 1; % state-specific fixed effects
[ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,x,N,T,model);

result2 = lmarginal_static_panel(ywith,xwith,W,N,T,info); 

fprintf(1,'region fixed effects: log marginal likelihoods and model probabilities \n');
in.cnames = strvcat('log-marginal','model probs');
in.rnames = strvcat('model','slx','sdm','sdem');
in.width = 10000;
in.fmt = '%10.4f';
out2 = [result2.lmarginal result2.probs];
mprint(out2,in);

model = 3; % state- and time-specific effects
[ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,x,N,T,model);

result3 = lmarginal_static_panel(ywith,xwith,W,N,T,info); 

fprintf(1,'region and time fixed effects: log marginal likelihoods and model probabilities \n');
in.cnames = strvcat('log-marginal','model probs');
in.rnames = strvcat('model','slx','sdm','sdem');
in.width = 10000;
in.fmt = '%10.4f';
out3 = [result3.lmarginal result3.probs];
mprint(out3,in);

% =============================================
% estimate SDM and SDEM models

ndraw = 3000;
nomit = 1000;

prior.model = 3;
prior.novi_flag = 1;
result1 = sdm_panel_FE_g(y,x,W,T,ndraw,nomit,prior);
vnames = strvcat('cclaims','jobposts','athome');
prt_panel(result1,vnames);

result2 = sdem_panel_FE_g(y,x,W,T,ndraw,nomit,prior);
vnames = strvcat('cclaims','jobposts','athome');
prt_panel(result2,vnames);