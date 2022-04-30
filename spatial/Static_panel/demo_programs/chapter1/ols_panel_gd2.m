% file ols_panel_gd2.m demo file2
clear all;
n = 200;
t = 10;

k = 2;
x = randn(n*t,k);
beta = ones(k,1);
sige = 0.1;
evec = randn(n*t,1)*sqrt(sige);

tts = (1:n)*(1/n);
SFE = kron(ones(t,1),tts');
ttt = (1:t)*(1/t);
TFE = kron(ttt',ones(n,1));

y = (x*beta + SFE + TFE + evec);

ndraw = 2500;
nomit = 500;
prior.novi_flag = 1; % homoscedastic model v_{it} = 1
prior.model = 3;     % model with fixed effects for regions and time period
prior.fe = 1;
result1 = ols_panel_FE_g(y,x,t,ndraw,nomit,prior);
vnames = strvcat('y','x1','x2');
prt_panel(result1,vnames);

tt = 1:n;
subplot(2,1,1),
plot(tt,0.5*result1.con+result1.sfe,'o',tt,tts,'+');
xlabel('N \times T observations');
ylabel('region-specific FE');
legend('estimated','true');
subplot(2,1,2),
tt=1:t;
plot(tt,0.5*result1.con+result1.tfe,'o',tt,ttt,'+');
xlabel('N \times T observations');
ylabel('time-specific FE');
legend('estimated','true');