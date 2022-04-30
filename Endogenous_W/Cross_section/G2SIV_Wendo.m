function results=G2SIV_Wendo(y,W,Zn,x1,x2,n)

%I. ESTIMATION of COEFFICIENTS

Gamma=x2\Zn; 
ep_hat=Zn-x2*Gamma;
Hn=[W*y x1 ep_hat];
Xn=[x1 x2(:,2:end)]; %constante always in the first column
Xn_tild=[x1(:,2:end) x2(:,2:end)]; 
Pn=x2*inv(x2'*x2)*x2';



%A. TFind initial values 
Qn_temp=[W*Xn W*Zn Xn_tild Zn]; % Use of Xn_tild insterad of Xn to avoid a multicollineary problem
kappa_2SIV=inv(Hn'*Qn_temp*inv(Qn_temp'*Qn_temp)*Qn_temp'*Hn)*Hn'*Qn_temp*inv(Qn_temp'*Qn_temp)*Qn_temp'*y;

%kappa_SIV=[rho beta delata]
rho_ini=kappa_2SIV(1,1);
beta_ini=kappa_2SIV(2:size(x1,2)+1,1);
delta_ini=kappa_2SIV(size(x1,2)+2:end,1);




%Constuct instruments from initial values
Gn_ini=W*inv(eye(n)-rho_ini*W);
Qn_ini=[Gn_ini*Xn Gn_ini*Zn Xn_tild Zn];

%Approximation  of PIn
sig2_ep=inv(n)*ep_hat'*ep_hat;
eta_ini=y-rho_ini*W*y-x1*beta_ini-ep_hat*delta_ini;
sig2_eta_ini=inv(n)*eta_ini'*eta_ini;
PIn_ini=sig2_eta_ini*eye(n)+delta_ini'*sig2_ep*delta_ini*Pn;



%G2SIV estimation from intial values
A1=inv(Hn'*inv(PIn_ini)*Qn_ini*inv(Qn_ini'*inv(PIn_ini)*Qn_ini)*Qn_ini'*inv(PIn_ini)*Hn);
kappa_G2SIV_ini=A1*Hn'*inv(PIn_ini)*Qn_ini*inv(Qn_ini'*inv(PIn_ini)*Qn_ini)*Qn_ini'*inv(PIn_ini)*y;

%B. To avoid dependence wrt initial values ,we optimize a second time using
%as initial values the values obtained as solution of the previous step.

%kappa_SIV=[rho beta delata]
rho_ini2=kappa_G2SIV_ini(1,1);
beta_ini2=kappa_G2SIV_ini(2:size(x1,2)+1,1);
delta_ini2=kappa_G2SIV_ini(size(x1,2)+2:end,1);

% %Construction of Instruments
Gn_ini2=W*inv(eye(n)-rho_ini2*W);
Qn_ini2=[Gn_ini2*Xn Gn_ini2*Zn Xn_tild Zn];

%Approximation of PIn
eta_ini2=y-rho_ini2*W*y-x1*beta_ini2-ep_hat*delta_ini2;
sig2_eta_ini2=inv(n)*eta_ini2'*eta_ini2;
Pn=x2*inv(x2'*x2)*x2';
PIn_ini2=sig2_eta_ini2*eye(n)+delta_ini2'*sig2_ep*delta_ini2*Pn;

%G2SIV estimation
kappa_G2SIV=inv(Hn'*inv(PIn_ini2)*Qn_ini2*inv(Qn_ini2'*inv(PIn_ini2)*Qn_ini2)*Qn_ini2'*inv(PIn_ini2)*Hn)*Hn'*inv(PIn_ini2)*Qn_ini2*inv(Qn_ini2'*inv(PIn_ini2)*Qn_ini2)*Qn_ini2'*inv(PIn_ini2)*y;


%C. Collect the final estimated values of the parameters 
rho=kappa_G2SIV(1,1);
beta=kappa_G2SIV(2:size(x1,2)+1,1);
delta=kappa_G2SIV(size(x1,2)+2:end,1);


Gn=W*inv(eye(n)-rho*W);
Qn=[Gn*Xn Gn*Zn Xn_tild Zn];

%Approximation of PIn
eta=y-rho*W*y-x1*beta-ep_hat*delta;
sig2_eta=inv(n)*eta'*eta;
PIn=sig2_eta*eye(n)+delta'*sig2_ep*delta*Pn;




sig_vep=sig2_ep*delta;
sig2_v=sig2_eta_ini+delta'*sig2_ep*delta;
coeffcorr=sig_vep*inv(sqrt(sig2_ep)*sqrt(sig2_v));

theta=[rho;beta;Gamma;sig2_ep;sig2_v;coeffcorr];


%II. Estimation of Var cov matrix  of kappa_G2SIV

Un=[Gn*(x1*beta+ep_hat*delta) x1 ep_hat];
sig_GIV=inv(n)*inv(Un'*inv(PIn)*Qn*inv(Qn'*inv(PIn)*Qn)*Qn'*inv(PIn)*Un);

results.sig_GIV=sig_GIV;
results.theta=theta;
results.kappa_G2SIV=kappa_G2SIV;























