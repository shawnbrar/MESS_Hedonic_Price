function results=gmm_sar_hetero_aug_p2_mc(y,x,W,Wy,Q,P,para_ini,We,sigma,h)
options=optimset('fminunc');
% options=optimset(options, 'MaxFunEvals',2000, 'LargeScale','off','TolFun', 0.00001, 'TolX',0.00001,'Display','off');
options=optimset(options, 'MaxFunEvals',2000, 'LargeScale','off','TolFun', 0.0000001, 'TolX',0.00001,'Display','off'); 
f_sar1gmm_oc=@(r) f_sar1gmm(r,y,x,Wy,Q,P,We);

[para]=fminunc(f_sar1gmm_oc,para_ini,options);
lam=para(1);
bet=para(2:end);
n=length(y);
In=eye(n);
A=In-lam*W;
%%%%%%%%%%%%%%%%%
G=W/A;


Q=[G*x(:,1:end-1)*bet(1:end-1,1), x(:,1:end-1),h]; %The last element of x (the predictor) is endogenous. 
P=G-diag(diag(G));

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

f_sar1gmm_oc=@(r) f_sar1gmm(r,y,x,Wy,Q,P,We);

para1=fminunc(f_sar1gmm_oc,para,options);
lam=para1(1);
bet=para1(2:end);
A=In-lam*W;
%%%%%%%%%%%%%%%%%
G=W/A;

Q=[G*x(:,1:end-1)*bet(1:end-1,1), x(:,1:end-1),h];
P=G-diag(diag(G));

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

% Computation of the variance-covariance matrix
D=zeros(size(Q,2)+1,size(x,2)+1); %k+1 parameters and q+1 moments
Ps=P+P';
SPsG=sigma*Ps*G;
D(1,1)=trace(SPsG);
D(2:end,1)=(G*x*bet)'*Q;
D(2:end,2:end)=Q'*x;
results.par=para1;
q=size(x,2)+1;
results.var=eye(q)/(D'*We*D);
results.resid=y-lam*Wy-x*bet;
% Attention, D in our notation corresponds to D' in Lee JoE 2007.