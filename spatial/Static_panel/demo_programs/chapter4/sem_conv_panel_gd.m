% sem_conv_panel_g demo file
clear all;
rng(19203040);
[unclaims,b] = xlsread('../demo_data/weekly.xlsx',1);
% read data from sheet 1 of Excel spreadsheet
% growth rate of unemployment claims 2019-2020 from same week, previous year
snames = strvcat(b(2:end,1)); % 48 state names
tnames = strvcat(b(1,2:end)); % 51 week labels
[N,T] = size(unclaims);
[jobposts,b] = xlsread('../demo_data/weekly.xlsx',2);
% read data from sheet 2 of Excel spreadsheet
% change in job posts from 1st week of 2020
[athome,b] = xlsread('../demo_data/weekly.xlsx',3);
% read data from sheet 3 of Excel spreadsheet
% growth rate of percent population at home
% 2019-2020 from same week, previous year
[a,b] = xlsread('../demo_data/Wcont48.xlsx');
% 48 x 48 contiguity matrix for states
Wcontiguity = normw(a);
% state-to-state commodity flows, 2017
[a,b] = xlsread('../demo_data/cflows_2017.xlsx');
% set main diagonal (intrastate flows) to zero
diaga = diag(a);
W = a - diag(diaga);
Wcom_flows = normw(W); % row-normalize
% eliminate small elements
for i=1:N
    for j=1:N
        if Wcom_flows(i,j) < 0.005
            Wcom_flows(i,j) = 0;
        end
    end
end

Wcom_flows = normw(Wcom_flows);

y = vec(unclaims);
x = [vec(jobposts) vec(athome)];

vnames = strvcat('y=unclaims','jobposts','athome');

Wmatrices = [Wcontiguity Wcom_flows];

ndraw = 25000;
nomit = 5000;
prior.model = 3;
prior.thin = 4;
result1 = sem_conv_panel_g(y,x,Wmatrices,N,T,ndraw,nomit,prior);
prt_panel(result1,vnames);

result2 = sdem_conv_panel_g(y,x,Wmatrices,N,T,ndraw,nomit,prior);
prt_panel(result2,vnames);




