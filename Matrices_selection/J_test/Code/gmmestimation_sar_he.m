function results=gmmestimation_sar_he(y,x,W,Q,P1)
%ATTENTION: this fonction only consider two step GMM estimation
% para is the initial set of parameters [\lambda \beta]
% P is the initial matrix used for the quadratic moment
% W is the interaction matrix
% n is the number of observations
% k is the number of columns in x
% Q is the matrix of instruments except the quadratic one (in the code,
% only one quadratic moment is used)
% We is the initial weight matrix used for the GMM approach (g'Wg)
% nlag=round(size(Y,1)^(1/3));
options = optimset;
% Modify options setting
options = optimset(options,'Display', 'off');
options = optimset(options,'FunValCheck', 'on');
%options = optimset(options,'LargeScale', 'off');
%options = optimset(options,'Diagnostics', 'off');

k=size(x,2);
para0=[0; zeros(k,1)];
% options =optimset(options, 'GradObj','on');
[n k]=size(x);
In=eye(n);
%  I used the mean of the g() function in the minimization, as you
%  suggested.
Wy=W*y;
q=size(Q,2);
We=eye(q+1);
% [para, fval, exitflag, output]=fmincon(@f_gmm_sar,para0,[],[],[],[],lb,ub,[],options,y,x,W,Q,P,We);
f_gmm_sar_he_oc=@(r) f_gmm_sar_he(r,y,x,W,Q,P1,We);
[para]=fminunc(f_gmm_sar_he_oc,para0,options);
% f0=feval(moment,paraest,3,y,x,W,Q,We);
% q=size(f0,1);
resid= (eye(n)-para(1)*W)*y-(x*para(2:end));

s=resid.*resid;
Sigm=diag(s);
%  Use of the best P matrix, namely P=G-diag(G), with G = W (I-\lambda
%  W)^{-1}
G=W*eye(n)/(eye(n)-para(1)*W);
P=G-diag(diag(G));
% Use of the best Q matrix, namely [GXb, X]
Q1=[G*x*para(2:end), x];
K=size(Q1,2);
Sp1=Sigm*P;
Sp1s=Sigm*(P+P'); 
% Attention, there is a typo in the Lin and Lee paper. The right formula is
% the one implemented here : tr(SPS(P^s)) and not tr(SP(SP)^s)
SpSp=Sp1.*Sp1s';
tr=sum(sum(SpSp));
qpsq=Q1'*Sigm*Q1;
Omeghe=zeros(K+1);
Omeghe(1,1)=tr;
Omeghe(2:end,2:end)=qpsq;
W1h=eye(K+1)/Omeghe;
count=1;
para_old=para0;
while sum(abs(para-para_old))> 0.00001 && count <20
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
Sigm=diag(res.*res);
Sp1=Sigm*P;
Sp1s=Sp1+Sp1';
% Sp2=sigma*p2;
% Sp2s=Sp2+Sp2';
delta11=sum(sum(Sp1.*Sp1s',2),1);
% delta21=sum(sum(Sp2.*Sp1s',2),1);
% delta22=sum(sum(Sp2.*Sp2s',2),1);
% delta12=sum(sum(Sp1.*Sp2s',2),1);

delta2= Q'*Sigm*Q;
% omega=zeros(2+size(Q1,2));
% omega(1:2,1:2)=[delta11 delta12; delta21 delta22];
% omega(3:end,3:end)=delta2;
omega=zeros(1+size(Q,2));
omega(1,1)=delta11;
omega(2:end,2:end)=delta2;
We=eye(1+size(Q,2))/(omega);
count=count+1;
end

f_gmm_sar_he_oc=@(r) f_gmm_sar_he(r,y,x,W,Q1,P1,W1h);
[para1]=fminunc(f_gmm_sar_he_oc,para,options);

% [para2,fval,exitflag,output,lambda,grad,hessian]=fmincon(@f_gmm_sar,para,[],[],[],[],lb,ub,[],options,y,x,W,Q1,P1,W1h);

resid1=(eye(n)-para1(1)*W)*y-(x*para1(2:end));

s=resid1.*resid1;
Sigm=diag(s);
%  Use of the best P matrix, namely P=G-diag(G), with G = W (I-\lambda
%  W)^{-1}
G=W*eye(n)/(eye(n)-para1(1)*W);
P1=G-diag(diag(G));
% Use of the best Q matrix, namely [GXb, X]
Q1=[G*x*para1(2:end), x];
K=size(Q1,2);
Sp1=Sigm*P1;
Sp1s=Sigm*(P1+P1'); 
SpSp=Sp1.*Sp1s';
tr=sum(sum(SpSp));
qpsq=Q1'*Sigm*Q1;
Omeghe=zeros(K+1);
Omeghe(1,1)=tr;
Omeghe(2:end,2:end)=qpsq;
W1h=eye(K+1)/Omeghe;

results.rho=para1(1);
 results.lambdah=para(1);
results.beta=para1(2:end);
 results.betah =para1(2:end);

% Computation of the variance-covariance matrix

results.resid=resid1;


G1h=W/(eye(n)-para1(1)*W);
Q1h=[G1h*x*para1(2:end) x];
Qpx=Q1h'*x;
P1h=G1h-diag(diag(G1h));
P2h=P1h+P1h';
Sp2=Sigm*P2h;
    trh=trace(Sp2*G1h);
    % omputation of D in Lee Lin paper
    Dh=zeros(K+1,k+1);
    Dh(1,1)=trh;
    Dh(2:end,1)=Q1h'*G1h*x*para1(2:end);
    Dh(2:end,2:end)=Qpx;
    
    results.Sigm=Sigm;
    xpxi = Dh'*W1h*Dh;
    cov=eye(k+1)/xpxi;
    results.cov=cov;
%     lik=-n/2*log(2*pi)-log(det(Sigm1))/2 + log(det(eye(n)-para2(1)*W))-residh'*inv(Sigm1)*residh;
%     results.lik=lik;
    tmph=diag(cov);
    results.tstat=[results.lambdah;results.betah]./sqrt(tmph);
    results.y=y;
    results.yhat=((eye(n)-para1(1)*W))\(x*para1(2:end));
    results.desc_tstat='first tstat is for \lambda then all following are for \beta';
    results.pvalues='not reported since we assume heteroskedastic and non normal error terms';
 
end
               
            

