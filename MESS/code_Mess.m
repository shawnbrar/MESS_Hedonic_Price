clear 
% load FDI_09
load FDI_09_1
load market_pot
mark = market_pot_debarsy;
% market potential is created as the log GDP of 183 countries weighted by
% the inverse distance (see Blonigen et al (2007) for further details.
%  It's based on the log GDP expressed in US of 2000

ly_con=log(FDI_09_con_Stock);% Log of dependent variable in constant USD of 2000
lgdp_con=log(GDP_cons_09); % log of GDPexpressed in constant USD of 2000
lpop=log(pop); %Log of population
ldistcap=log(distcap*1000);
mrp=tariff_rate_mean_all_product; %Some measure of tariff rates 
wmrp=tariff_rate_wmean_all_product;% Some other measure of tariff rates
lat=Lat; %Latitute of countries
n=size(lat,1);
long=Long;% Longitude
oecd=OECDDummy; % Oecd dummy
W=dinvarc(Lat,Long,1); % Create inverse distance matrix 
d=eigs(W);

Ws=W/d(1,1);

x_con=[ones(n,1) lgdp_con lpop oecd ldistcap wmrp mark];

info.lflag=0;
tot=strvcat('log FDI Constant', 'Constant', 'lgdp','lpop', 'oecd','distcap', 'Tariffs', 'market potential');
 res0=hwhite(ly_con,x_con);
prt(res0,tot);
% MESS estimation 
res1b=mess_nls(ly_con, x_con,Ws); % beta, alpha
 Q= [Ws*x_con(:,[2:end]) x_con];
 P=Ws;
    We=eye(size(Q,2)+1);
    resgm=gmmestimation_mess(ly_con,x_con,Ws,Q,P,We); % alpha, beta sigma2
    
    % % % % Computation of impacts For Mess  model by ML 
% % GDP
impact_mess_2 = expm(-res1b.alpha*Ws)*res1b.beta(2);
d2m = diag(impact_mess_2);
% % pop
impact_mess_3 = expm(-res1b.alpha*Ws)*res1b.beta(3);
d3m = diag(impact_mess_3);
% % OECD
impact_mess_4 = expm(-res1b.alpha*Ws)*res1b.beta(4);
d4m = diag(impact_mess_4);
% % dist
impact_mess_5 = expm(-res1b.alpha*Ws)*res1b.beta(5);
d5m = diag(impact_mess_5);
% %Tariff
impact_mess_6 = expm(-res1b.alpha*Ws)*res1b.beta(6);
d6m = diag(impact_mess_6);

% %MP
impact_mess_7 = expm(-res1b.alpha*Ws)*res1b.beta(7);
d7m = diag(impact_mess_7);

% Computation of ADE
ex=expm(-res1b.alpha*Ws);
% Computation of the numerator
% For GDP
num2 = 1/n*trace(impact_mess_2);

% For pop 
num3 = 1/n*trace(impact_mess_3);
% For OECD
num4 = 1/n*trace(impact_mess_4);
% For dis
num5 = 1/n*trace(impact_mess_5);
% For tariffs
num6 = 1/n*trace(impact_mess_6);
% % For MP
num7 = 1/n*trace(impact_mess_7);

res1=res1b; %USe QML results
% For GDP 
A2=[1/n*trace(-Ws*ex)*res1.beta(2) 1/n*trace(ex)];
S2 = [res1.covr(end-1,end-1) res1.covr(end-1,2)
    res1.covr(2,end-1) res1.covr(2,2)];
den2=(A2*S2*A2')^(0.5);
tstat2 = num2/den2;
prob2=norm_prb(tstat2);
% For pop 
A3=[1/n*trace(-Ws*ex)*res1.beta(3) 1/n*trace(ex)];
S3 = [res1.covr(end-1,end-1) res1.covr(end-1,3)
    res1.covr(3,end-1) res1.covr(3,3)];
den3=(A3*S3*A3')^(0.5);
tstat3 = num3/den3;
prob3=norm_prb(tstat3);
% For oecd 
A4=[1/n*trace(-Ws*ex)*res1.beta(4) 1/n*trace(ex)];
S4 = [res1.covr(end-1,end-1) res1.covr(end-1,4)
    res1.covr(4,end-1) res1.covr(4,4)];
den4=(A4*S4*A4')^(0.5);
tstat4 = num4/den4;
prob4=norm_prb(tstat4);
% For dis
A5=[1/n*trace(-Ws*ex)*res1.beta(5) 1/n*trace(ex)];
S5 = [res1.covr(end-1,end-1) res1.covr(end-1,5)
    res1.covr(5,end-1) res1.covr(5,5)];
den5=(A5*S5*A5')^(0.5);
tstat5 = num5/den5;
prob5=norm_prb(tstat5);
% For tarifs
A6=[1/n*trace(-Ws*ex)*res1.beta(6) 1/n*trace(ex)];
S6 = [res1.covr(end-1,end-1) res1.covr(end-1,6)
    res1.covr(6,end-1) res1.covr(6,6)];
den6=(A6*S6*A6')^(0.5);
tstat6 = num6/den6;
prob6=norm_prb(tstat6);
% For MP
A7=[1/n*trace(-Ws*ex)*res1.beta(7) 1/n*trace(ex)];
S7 = [res1b.covr(end-1,end-1) res1b.covr(end-1,7)
    res1b.covr(7,end-1) res1b.covr(7,7)];
den7=(A7*S7*A7')^(0.5);
tstat7 = num7/den7;
prob7=norm_prb(tstat7);

display('Results of average direct effects for MESS hetero by QML')
% se_avgdir=[num2 num3 num4 num5 num6 num7
%     den2 den3 den4 den5 den6 den7
%     prob2 prob3 prob4 prob5 prob6 prob7  ]
    pavgdir=[num2 num3 num4 num5 num6 num7
%     den2 den3 den4 den5 den6 
    prob2 prob3 prob4 prob5 prob6 prob7  ]

% % Let's now analyze the significance of some effects
% for GDP, indirect between Austria and Slovakia
%  This is the highest negative value since the 2 countries are really
%  close
eni=zeros(n,1);
enj=zeros(n,1);
eni(31,:)=1;
enj(2,:)=1;
% Computation of the numerator of the statistic
num = eni'*expm(-res1b.alpha*Ws)*enj*res1b.beta(2);
An=[(-eni'*Ws*expm(-res1b.alpha*Ws)*enj)*res1b.beta(2) eni'*expm(-res1b.alpha*Ws)*enj ];
cov =res1b.cov; % load of the homoskedastic consistent cov matrix
Bn=[res1b.cov(end-1,end-1) res1b.cov(end-1,2)
    res1b.cov(2,end-1) res1b.cov(2,2)]; % the spatial autoregressive parameter is in the last row (column)
VC=An*Bn*An';
SE = sqrt(VC);
imp_gdp_svk_aut=num/SE;
display('Prob that AUT and SVK are significantly different from 0 for GDP')
prob_gdp=norm_prb(imp_gdp_svk_aut)


%%% SAR estimation and inference on impacts
info.lflag=0;
%%%%%%%%%%%%%%%
%ATTENTION%%%
%%%%%%%%%

%open the function sar: open sar
% On line 206 add results.cov=hessi(1:nvar+1, 1:nvar+1);

res_sar=sar(ly_con, x_con,Ws,info);
prt(res_sar,tot)
In=eye(n);
% % Impact GDP
Si=In/(In-res_sar.rho*Ws);
impact_gdp = Si*res_sar.beta(2);
d2m = diag(impact_gdp);
fprintf('ADE for GDP is %f2.3\n',mean(d2m))

i2r=sum(Si-diag(diag(Si)),2)*res_sar.beta(2);
i2c=sum(Si-diag(diag(Si)),1)*res_sar.beta(2);

fprintf('Spillin effect for GDP is %f2.3\n',mean(i2r));
fprintf('Spillout effect for GDP is %f2.3\n',mean(i2r));

% Inference on ADE for GDP

=res_sar.cov;
Deriv=zeros(size(x_con,2)+1,1);
Deriv(2,1)=mean(diag(Si));
Deriv(end,1)=mean(diag(Si*Ws*Si))*res_sar.beta(2);
num=mean(d2m);
den=(Deriv'*cov*Deriv)^0.5;

t_stat_gdp=num/den;
prob=norm_prb(t_stat_gdp);

%% Inference on ATE for GDP 
