% sem_conv_panel_g demo file
% using the function to produce a set of cross-sectional model estimates
clear all;
rng(19203040);
[unclaims,b] = xlsread('../demo_data/weekly.xlsx',1);
% read data from sheet 1 of Excel spreadsheet
% growth rate of unemployment claims 2019-2020 from same week, previous year
snames = strvcat(b(2:end,1)); % 48 state names
tnames = strvcat(b(1,2:end)); % 51 week labels
[N,T] = size(unclaims);
T = 1; % only use the cross-section from the first time period
unclaims = unclaims(:,T);
[jobposts,b] = xlsread('../demo_data/weekly.xlsx',2);
% read data from sheet 2 of Excel spreadsheet
% change in job posts from 1st week of 2020
jobposts = jobposts(:,T);
[athome,b] = xlsread('../demo_data/weekly.xlsx',3);
% read data from sheet 3 of Excel spreadsheet
% growth rate of percent population at home
% 2019-2020 from same week, previous year
athome = athome(:,T);
[a,b] = xlsread('../demo_data/Wcont48.xlsx');
% 48 x 48 contiguity matrix for states
W = normw(a);

y = vec(unclaims);
x = [ones(N,1) vec(jobposts) vec(athome)];

vnames = strvcat('y=unclaims','constant','jobposts','athome');

ndraw = 2500;
nomit = 500;
prior.model = 0;
result1 = sem_panel_FE_g(y,x,W,T,ndraw,nomit,prior);
prt_panel(result1,vnames);





