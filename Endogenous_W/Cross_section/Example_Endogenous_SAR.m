% PURPOSE: An example of using the spatial panel estimation
% and spatial autocorrelation tests toolbox        
%Here, Wn is endogenous
%---------------------------------------------------

clc
clear all;





%B. Loading the spatial weight matrix
load mat81q;
Wbis=mat81q;
Wd=(Wbis>0);% binary matrix (1 if connected and 0 otherwize)
clear mat81q;




N=size(Wd,1);



%II.Construction  SAR
nvar=1;
randn('state',14789874);    
x1=[ones(N,1) randn(N,nvar)]; %Always put the intercept in first column
b=ones(nvar+1,1); %Coeff variable x
rho=0.4; %autoccorelation mean part



%III. Endogeneous part
%A. Construction of error term V_n and epsilon_n
mu = [0 0];
sigma_ep_v=0.8;
sig2_v=1;
p2=1; % Dimension of Z
sig2_ep=1;

coeffcorr=sigma_ep_v/(sqrt(sig2_ep)*sqrt(sig2_v));
Sigma = [sig2_v coeffcorr; coeffcorr sig2_ep];
Gamma=[1 0.8]';

% One constructs x2 (IV of first step)
% x2 is used to construct Z, used to construct W.

randn('state',156488);    
x2=[ones(N,1) randn(N,1)]; 






%hwait = waitbar(0,'Sampling...');
incr=2;
% for incr=1:nsim
%    try
%waitbar(incr/nsim);
randn('state',46154185+incr);    
r = mvnrnd(mu,Sigma,N);
Vn=r(:,1);
Ep_i=r(:,2:end);
%B.Construction de Zn
Zn=kron(ones(1,p2),x2*Gamma)+Ep_i;
We=zeros(N);
for i=1:N
    for j=i+1:N
   We(i,j)=1/abs(Zn(i)-Zn(j));
    end
end
for j=1:N
    for i=j+1:N
   We(i,j)=1/abs(Zn(j)-Zn(i));
    end
end


%III. Construction of weight matrix
W=Wd.*We;
%W=We;
%Row-normalization
Wn = normw(W);   
%  Normalization KP (lemme1)
d=eig(W);
d1=max(d);
Wn=W/d1; 


%% Construction of data
In=eye(N);



A=In-rho*Wn;

Ai=In/A;


%Generate the DGP
y=Ai*(x1*b+Vn);



%% Regression


theta_0=[rho;b;Gamma;sig2_ep;sig2_v;coeffcorr];
vnames = strvcat('y','const','x1');

info.lflag=0;
info.ndraw = 2500;

%G2SIV, endogenous W
resultG2SIV_Wendo=G2SIV_Wendo(y,Wn,Zn,x1,x2,N);
[theta_0 resultG2SIV_Wendo.theta]
% pkg load optim 
% pk g load statistics

%QMLE, endogenous W
resultQMLE_Wendo=sar_Wendo(y,x1,Wn,Zn,x2,info);
%resultQMLE_Wendo = sar_Wendo(y,X,Wn,Zn,x2,info);
theta_hat=[resultQMLE_Wendo.rho;resultQMLE_Wendo.beta;resultQMLE_Wendo.Gamma;resultQMLE_Wendo.sig2_ep;resultQMLE_Wendo.sig2_v;resultQMLE_Wendo.coeffcorr];
%Comparison with the true value
tstat=resultQMLE_Wendo.Allparam./sqrt(diag(resultQMLE_Wendo.Allvarcov)); %[rho beta gamma sig2_ep sig2_v coeffcorr]
[theta_0 theta_hat tstat]
prt_Wendo(resultQMLE_Wendo,vnames);
pause
%QMLE, exogenous W
result = sar(y,x1,Wn,info); % maximum likelihood estimates
prt(result,vnames);
pause
%GMM, exogenous W
result3 = sar_gmm(y,x1,Wn,info);
prt(result3,vnames);


%catch ee 
%end
%end
%close(hwait);





