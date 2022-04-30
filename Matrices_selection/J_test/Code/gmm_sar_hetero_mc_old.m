function results=gmm_sar_hetero_mc_old(y,x,W,Wy,Q,P,para_ini,We)
% This function uses the omega matrix of Lin and Lee(2010) 
options=optimset('fminsearch');
% options=optimset(options, 'MaxFunEvals',2000, 'LargeScale','off','TolFun', 0.00001, 'TolX',0.00001,'Display','off');
options=optimset(options, 'MaxFunEvals',2000, 'TolFun', 0.0001, 'TolX',0.0001,'Display','off'); 

[para]=fminsearch('f_sar1gmm',para_ini,options,y,x,Wy,Q,P,We);
lam=para(1);
bet=para(2:end);
k=size(para,1);
n=length(y);
In=eye(n);
A=In-lam*W;
%%%%%%%%%%%%%%%%%
G=W/A;

Q=[G*x*bet, x];
P=G-diag(diag(G));
res=y-lam*Wy-x*bet;
sigma = diag(res.*res);
Sp1=sigma*P;
Sp1s=Sp1+Sp1';
% Sp2=sigma*p2;
% Sp2s=Sp2+Sp2';
delta11=sum(sum(Sp1.*Sp1s',2),1);
% delta21=sum(sum(Sp2.*Sp1s',2),1);
% delta22=sum(sum(Sp2.*Sp2s',2),1);
% delta12=sum(sum(Sp1.*Sp2s',2),1);

delta2= Q'*sigma*Q;
% omega=zeros(2+size(Q1,2));
% omega(1:2,1:2)=[delta11 delta12; delta21 delta22];
% omega(3:end,3:end)=delta2;
omega=zeros(1+size(Q,2));
omega(1,1)=delta11;
omega(2:end,2:end)=delta2;
We=eye(1+size(Q,2))/(omega);
%%%%%%%%%%%%%%%


count=1;
para_old=para_ini;
while sum(abs(para-para_old))> 0.0001 && count <20
    para_old=para;
para=fminsearch('f_sar1gmm',para_old,options,y,x,Wy,Q,P,We);
lam=para(1);
bet=para(2:end);
results.par=para;
A=In-lam*W;
%%%%%%%%%%%%%%%%%
G=W/A;

Q=[G*x*bet, x];
P=G-diag(diag(G));
res=A*y-x*bet;
sigma=diag(res.*res);
Sp1=sigma*P;
Sp1s=Sp1+Sp1';
% Sp2=sigma*p2;
% Sp2s=Sp2+Sp2';
delta11=sum(sum(Sp1.*Sp1s',2),1);
% delta21=sum(sum(Sp2.*Sp1s',2),1);
% delta22=sum(sum(Sp2.*Sp2s',2),1);
% delta12=sum(sum(Sp1.*Sp2s',2),1);

delta2= Q'*sigma*Q;
% omega=zeros(2+size(Q1,2));
% omega(1:2,1:2)=[delta11 delta12; delta21 delta22];
% omega(3:end,3:end)=delta2;
omega=zeros(1+size(Q,2));
omega(1,1)=delta11;
omega(2:end,2:end)=delta2;
We=eye(1+size(Q,2))/(omega);
count=count+1;
end


results.par=para;
results.resid=y-lam*Wy-x*bet;
results.sigma=sigma;
