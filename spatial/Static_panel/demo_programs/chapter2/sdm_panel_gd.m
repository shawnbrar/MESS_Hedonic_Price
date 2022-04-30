% sdm_panel_g demo file
clear all;
[uclaims,b] = xlsread('../demo_data/weekly.xlsx',1);
% read data from sheet 1 of Excel spreadsheet
% growth rate of unemployment 2019-2020 from same week, previous year
snames = strvcat(b(2:end,1)); % 48 state names
tnames = strvcat(b(1,2:end)); % 51 week labels
[N,T] = size(uclaims);
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

y = vec(uclaims);
x = [vec(jobposts) vec(athome)];

info.model = 3;
result1 = sdm_panel_FE(y,x,W,T,info);
vnames = strvcat('uclaims','jobposts','athome');
prt_panel(result1,vnames);

ndraw = 2500;
nomit = 500;
prior.novi_flag = 1;
prior.model = 3;
result2 = sdm_panel_FE_g(y,x,W,T,ndraw,nomit,prior);
prt_panel(result2,vnames);

prior2.rval = 5;
prior2.model = 3;
prior2.fe = 1;
result3 = sdm_panel_FE_g(y,x,W,T,ndraw,nomit,prior2);
prt_panel(result3,vnames,snames,tnames);

tt=1:T;
plot(tt,result3.tfe);
xlabel('weeks during 2020');
ylabel('time fixed effects estimates');

vmat = reshape(result3.vmean,N,T);
vmean = mean(vmat,1);
plot(tt,vmean,'o');
xlabel('Weeks during 2020');
ylabel('Mean of time v_{it} estimates');


