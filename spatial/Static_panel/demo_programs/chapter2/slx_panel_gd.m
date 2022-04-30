% slx_panel_gd 
clear all;
[unclaims,b] = xlsread('../demo_data/weekly.xlsx',1);
% read data from sheet 1 of Excel spreadsheet
% growth rate of unemployment 2019-2020 from same week, previous year
snames = strvcat(b(2:end,1)); % 48 state names
tnames = strvcat(b(1,2:end)); % 51 week labels
[N,T] = size(unclaims);
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

y = vec(unclaims);
x = [vec(jobposts) vec(athome)];
vnames = strvcat('y=unclaims','jobposts','athome');

info.model = 3;
ndraw = 2500;
nomit = 500;
info.novi_flag = 1;
result1 = slx_panel_FE_g(y,x,W,T,ndraw,nomit,info);
prt_panel(result1,vnames);

prior.novi_flag = 1;
prior.model = 3;
result2 = sdem_panel_FE_g(y,x,W,T,ndraw,nomit,prior);
prt_panel(result2,vnames);

betao = result1.bdraw(:,2);
betas = result2.bdraw(:,2);
beta_diff = betas - betao;

thetao = result1.bdraw(:,4);
thetas = result2.bdraw(:,4);
theta_diff = thetas - thetao;


[h1,f1,y1] = pltdens(betao+thetao);
[h2,f2,y2] = pltdens(betas+thetas);
[h3,f3,y3] = pltdens(beta_diff+theta_diff);

subplot(2,1,1),
plot(y1,f1,'.-r',y2,f2,'.-b');
ylabel('\beta + \theta posteriors');
xlabel('total effects values');
legend('\beta_o+\theta_o','\beta_s+\theta_s');
subplot(2,1,2),
plot(y3,f3,'.-g');
ylabel('Posterior for differences');
xlabel('total effects values');
zipi = find(y3 > 0);
line([0 0],[0 f3(zipi(1,1))]);
legend('total effects differences','zero');

% trapezoid rule integration
sum_all = trapz(y3,f3);
sum_positive = trapz(y3(zipi,1),f3(zipi,1));
prob = sum_positive/sum_all

