function results=gmm_sar_hetero_aug_p1_mc(y,x,W,Wy,Q,P,para_ini,We)
options=optimset('fminsearch');
% options=optimset(options, 'MaxFunEvals',2000, 'LargeScale','off','TolFun', 0.00001, 'TolX',0.00001,'Display','off');
options=optimset(options, 'MaxFunEvals',2000, 'FunValCheck', 'on','MaxIter', 1000,'TolFun', 0.0001, 'TolX',0.0001,'Display','off'); 
f_sar1gmm_oc=@(r) f_sar1gmm(r,y,x,Wy,Q,P,We);
[para]=fminsearch(f_sar1gmm_oc,para_ini,options);
lam=para(1);
bet=para(2:end);
n=length(y);
In=eye(n);
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
omega_old=zeros(1+size(Q,2));
omega_old(1,1)=delta11;
omega_old(2:end,2:end)=delta2;
if abs(det(omega_old))<1e-6
    omega_old=eye(1+size(Q,2));
end
We=eye(1+size(Q,2))/(omega_old);
%%%%%%%%%%%%%%%

count=1;
para_old=para_ini;
while sum(abs(para-para_old))> 0.0001 && count <50
    para_old=para;
f_sar1gmm_oc=@(r) f_sar1gmm(r,y,x,Wy,Q,P,We);
    
para=fminsearch(f_sar1gmm_oc,para_old,options);
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
if abs(det(omega))<1e-6 && abs(det(omega_old))<1e-6
    omega=eye(1+size(Q,2));
elseif abs(det(omega))<1e-6
    omega=omega_old;  
end

We=eye(1+size(Q,2))/(omega);
count=count+1;
end
% Computation of the variance-covariance matrix
D=zeros(size(Q,2)+1,size(x,2)+1); %k+1 parameters and q+1 moments
Ps=P+P';
SPsG=sigma*Ps*G;
D(1,1)=trace(SPsG);
D(2:end,1)=Q'*(G*x*bet);
D(2:end,2:end)=Q'*x;
results.par=para;
q=size(x,2)+1;
results.var=eye(q)/(D'*We*D);
results.std=sqrt(diag(results.var));
results.resid=y-lam*Wy-x*bet;
results.sigma=sigma;
results.count=count;
% Attention, D in our notation corresponds to D' in Lee JoE 2007.