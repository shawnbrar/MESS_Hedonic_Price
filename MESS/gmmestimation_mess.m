function results=gmmestimation_mess(y,x,W,Q,P,We)
% Q is the matrix of linear moments [Wx x]
% P is the matrix used for quadratic moments. 
% We is the initial weight matrix for the GMM approach. It's updated in the
% second step. 
[n k]=size(x);
% nlag=round(size(Y,1)^(1/3));
options=optimset('fminunc');
options=optimset(options, 'MaxFunEvals',2000);
options.LargeScale='off';
% options =optimset(options, 'GradObj','on');

%  I used the mean of the g() function in the minimization, as you
%  suggested.
res=ols(y,x);
para0=[0;zeros(size(res.beta,1),1)];
f_messgmm_oc=@(r) f_messgmm(r,y,x,W,Q,P,We);
[para]=fminunc(f_messgmm_oc,para0,options);
paraest=para;
% f0=feval(moment,paraest,3,y,x,W,Q,We);
% q=size(f0,1);
resid=expm(paraest(1)*W)*y-(x*paraest(2:end));
sige=resid'*resid/n;
sige2=sige^2;
%  definition of the size of the weighting matrix for gmm
K = size(Q,2)+1; % 1 quadratic moment
Omeg=zeros(K);
% Homoskedastic case
Ps=P+P';
vPs=vec(Ps);
vWs=vec(W+W');
tr=vPs'*vPs;
qpq=Q'*Q;
Omeg(1,1)=sige2/(2*n)*tr;
Omeg(2:end,2:end)=sige/n*qpq;

W1=eye(K)/Omeg;
f_messgmm_oc=@(r) f_messgmm(r,y,x,W,Q,P,W1);
[para1]=fminunc(f_messgmm_oc,paraest,options);

% f0=feval(moment,paraest,3,y,x,W,Q,We);
% q=size(f0,1);
residh=expm(para1(1)*W)*y-(x*para1(2:end));
sigeh=residh'*residh/n;
sige2h=sigeh^2;
Omeg=zeros(K);
% Homoskedastic case
Omeg(1,1)=sige2h/(2*n)*tr;
Omeg(2:end,2:end)=sige/n*qpq;

W1=eye(K)/Omeg;
%  Refinement of the homoskedastic case.
f_messgmm_oc=@(r) f_messgmm(r,y,x,W,Q,P,W1);
[para1ho]=fminunc(f_messgmm_oc,para1,options);
resid1=expm(para1ho(1)*W)*y-(x*para1ho(2:end));
sige1ho=resid1'*resid1/n;
Omeg=zeros(K);
Omeg(1,1)=sige2h/(2*n)*tr;
Omeg(2:end,2:end)=sige/n*qpq;


W1ho=eye(K)/Omeg;
%  Refinement of the homoskedastic case.
f_messgmm_oc=@(r) f_messgmm(r,y,x,W,Q,P,W1ho);
[para1ho]=fminunc(f_messgmm_oc,para1ho,options);
results.alpha=para1ho(1);
results.beta=para1ho(2:end);
resid1=expm(para1ho(1)*W)*y-(x*para1ho(2:end));
sige1ho=resid1'*resid1/n;
results.sige=sige1ho;
Omeg=zeros(K);
Omeg(1,1)=sige2h/(2*n)*tr;
Omeg(2:end,2:end)=sige/n*qpq;


W1ho=eye(K)/Omeg;


vW2=vPs'*vWs;
tr1=0.5*(vW2^2)/tr;
xb=x*results.beta;
wxb=W*xb;
qpwxb=Q'*wxb;
q1=size(Q,2);
% Computation of the COV matrix
qpqi=eye(q1)/qpq;
qpx=Q'*x;

G1=zeros(K,k+1);
G1(1,1)=sige1ho/2*vW2;
G1(2:end,1)=qpwxb;
G1(2:end,2:end)=-qpx;
delta = 1/n*G1'*W1ho*G1;
results.delta1=delta;
delta1i=eye(k+1)/delta;
temp1=diag(delta1i);
results.tstat_homo=[results.alpha;results.beta]./sqrt(temp1);

% Homoskedastic case and non normality
% It does not change anything since P has 0 diagonal. Thus, the matrix V is
% the same as for the normal case (see p.18 of the revised version).

% Heteroskedasticity robust version
s=resid1.*resid1; %Selection of the last estimated residuals since it's the most efficient.
Sigm=diag(s); % Heteroskedastic matrix.
Sd=sqrt(Sigm); % square root of the heteroskedastic matrix.
Omeghe=zeros(K);
Swpw=Sd*Ps*Sd;
omh=vec(Swpw);
qpsq=Q'*Sigm*Q;
Omeghe(1,1)=1/(2*n)*omh'*omh;
Omeghe(2:end,2:end)=1/n*qpsq;
W1h=eye(K)/Omeghe;
f_messgmm_oc=@(r) f_messgmm(r,y,x,W,Q,P,W1h);
[para2]=fminunc(f_messgmm_oc,para1ho,options);
% f0=feval(moment,paraest,3,y,x,W,Q,We);
% q=size(f0,1);
residhe=expm(para2(1)*W)*y-x*para2(2:end);
s=residhe.*residhe;
Sigmhe=diag(s);
Sd=sqrt(Sigmhe);
Omeghe=zeros(K);
omh=vec(Sd*Ps*Sd);
qpsq=Q'*Sigmhe*Q;
Omeghe(1,1)=0.5/n*omh'*omh;
Omeghe(2:end,2:end)=1/n*qpsq;
W1he=eye(K)/Omeghe;
% Refinement of the heteroskedastic case
f_messgmm_oc=@(r) f_messgmm(r,y,x,W,Q,P,W1he);
[para3]=fminunc(f_messgmm_oc,para2,options);
residhe1=expm(para3(1)*W)*y-(x*para3(2:end));
s=residhe1.*residhe1;
Sigmhe1=diag(s);
Sd=sqrt(Sigmhe1);
Omeghe=zeros(K);
omh=vec(Sd*Ps*Sd);
qpsq=Q'*Sigmhe1*Q;
Omeghe1=zeros(K);
Omeghe1(1,1)=1/(2*n)*omh'*omh;
Omeghe1(2:end,2:end)=1/n*qpsq;
W1he1=eye(K)/Omeghe1;

% Refinement of the heteroskedastic case
f_messgmm_oc=@(r) f_messgmm(r,y,x,W,Q,P,W1he1);
[para4]=fminunc(f_messgmm_oc,para3,options);
residhe1=expm(para4(1)*W)*y-(x*para4(2:end));
s=residhe1.*residhe1;
Sigmhe1=diag(s);
Sd=sqrt(Sigmhe1);
Omeghe=zeros(K);
omh=vec(Sd*Ps*Sd);
qpsq=Q'*Sigmhe1*Q;
Omeghe1=zeros(K);
Omeghe1(1,1)=1/(2*n)*omh'*omh;
Omeghe1(2:end,2:end)=1/n*qpsq;
W1he1=eye(K)/Omeghe1;
results.alphah=para4(1);
results.betah=para4(2:end);

% Computation of the variance-covariance matrix
xb=x*para4(2:end);
wxb=W*xb;
qpwxb=Q'*wxb;

Sii=inv(Sigmhe1);
SigmW=Sii*W;
SigmWs=SigmW + SigmW';
omss=vec(Sd*SigmWs*Sd);
G1=zeros(K,k+1);
G1(1,1)=0.5*omh'*omss;
G1(2:end,1)= qpwxb;
G1(2:end,2:end)=-qpx;
delta_het=1/n*G1'*W1he1*G1; %Premultiplication by 1/n to compensate for the n in the W1he1 matrix. 
% Like this, we look at the distribution of the parameters and not of
% sqrt(n) of the parameters.
deltai_het=eye(k+1)/delta_het;
temp_het=diag(deltai_het);
results.tstat_hetero=[results.alphah;results.betah]./sqrt(temp_het);
% 
% 
% delta_h=zeros(k+1);
% delta_h(1,1)=tr1h+qpwxb'*(qpqi*qpwxb);
% delta_h(2:end,1)=-(qpx)'*(qpqi*qpwxb);
% delta_h(1,2:end)=(delta_h(2:end,1))';
% delta_h(2:end,2:end)=qpx'*(qpqi*qpx);
% results.delta=delta_h;
% deltai=eye(k+1)/delta_h;
% tmph=diag(deltai);
% results.tstat_hetero=[results.alphah;results.betah]./sqrt(tmph);