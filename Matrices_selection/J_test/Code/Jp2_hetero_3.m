function [stat1, stat2,stat3,signi1]=Jp2_hetero_3(y,W_1,W_2,W_3,x,res_1,res_2,res_3,B)
% This code is for the first predictor.
% We don't need to instrument the predictor since it's exogenous. 
% So it's just an IV to instrument the Wy using the W which is assumed the
% true one
rho_w1=res_1.rho;
beta_w1=res_1.beta;
res_w1=res_1.resid;

rho_w2=res_2.rho;
beta_w2=res_2.beta;
res_w2=res_2.resid;

rho_w3=res_3.rho;
beta_w3=res_3.beta;
res_w3=res_3.resid;

% First need to compute the predictor for the tree models
[n,k]=size(x);
Con=[zeros(2,k+1), eye(2)]; %COnstraint that the predictors of the alternative models are not significant.
In=eye(n);
A1i=In/(In-rho_w1*W_1);

A2i=In/(In-rho_w2*W_2);
A3i=In/(In-rho_w3*W_3);


y1p=rho_w1*W_1*y+(x*beta_w1);
y2p=rho_w2*W_2*y+(x*beta_w2);
y3p=rho_w3*W_3*y+(x*beta_w3);

% Compute the statistic assuming that each model can be the true

% Computation of the sigma under the trhee possible null hypotheses
% We do not use that anymore because results are better when the variance
% is estimated from the augmented model
% sigma1=diag(res_w1.*res_w1);
% sigma2=diag(res_w2.*res_w2);
% sigma3=diag(res_w3.*res_w3);


% Explanatory matrix under the three different null

X1=[x, y2p y3p];% Null is that W_1 is the true
X2=[x, y1p y3p]; % Null is that W_2 is the true
X3=[x, y1p y2p]; % Null is that W_3 is the true

% Linear moments for estimation of the three augmented models 


if sum(sum(W_1,2),1)==n
W1x=W_1*x(:,2:end);
W2x=W_2*x(:,2:end);
W3x=W_3*x(:,2:end);
else
 W1x=W_1*x;
W2x=W_2*x;
W3x=W_3*x;
end
h1=[W_1*W1x W1x];
h2=[W_2*W2x W2x];
h3=[W_3*W3x W3x];

% Moments conditions for the estimation of null model
Q1=[h1 x];
Q2=[h2 x];
Q3=[h3 x];

% Moment conditions for the estimation of augmented model

Q1a=[h1 X1];
Q2a=[h2 X2];
Q3a=[h3 X3];



% Computation of the quadratic moment for the three augmented models
p1 = W_1; % Homoskedastic case. We use the quadratic moment of Lee, Journal of Econometrics, 2007.
p2 = W_2;
p3 = W_3;


W1y=W_1*y; 
W2y=W_2*y;
W3y=W_3*y;



We=eye(size(Q1a,2)+1); % One quadratic moment

% Computation of the augmented model for each null
para_ini=zeros(size(X1,2)+1,1); %k +1 parameters to estimate

results1=gmm_sar_hetero_aug_p1_mc(y,X1,W_1,W1y,Q1a,p1,para_ini,We);

results2=gmm_sar_hetero_aug_p1_mc(y,X2,W_2,W2y,Q2a,p2,para_ini,We);

results3=gmm_sar_hetero_aug_p1_mc(y,X3,W_3,W3y,Q3a,p3,para_ini,We);




 delt1=Con*[results1.par];
 delt2=Con*[results2.par];
 delt3=Con*[results3.par];



 varcov1=results1.var;
 varcov2=results2.var;
 varcov3=results3.var;


 
 
 Rv1Ri=eye(2)/(Con*varcov1*Con');
 Rv2Ri=eye(2)/(Con*varcov2*Con');
 Rv3Ri=eye(2)/(Con*varcov3*Con');



%  Computation of the three J test stats
 stat1=delt1'*Rv1Ri*delt1;
 stat2=delt2'*Rv2Ri*delt2;
 stat3=delt3'*Rv3Ri*delt3;

% Then compare the statistics 
[~,ind]=min([stat1, stat2, stat3],[],2)
J1_boot=zeros(B,3);
par_ini_n=zeros(size(x,2)+1,1);
Wen=eye(size(Q1,2)+1);
para_ini=zeros(size(X1,2)+1,1); %Initial value for lambda
% Then bootstrap the statistic and compare

if ind==1% The min values is for W1
for b=1:B
    res_w1_star=reechantillonage(res_w1);
    ystar=A1i*(x*beta_w1+res_w1_star);
   w1ys=W_1*ystar;
   w2ys=W_2*ystar;
   w3ys=W_3*ystar;
               RGMM2st=gmm_sar_hetero_mc(ystar,x, W_2,w2ys,Q2,p2,par_ini_n,Wen);
               RGMM3st=gmm_sar_hetero_mc(ystar,x, W_3,w3ys,Q3,p3,par_ini_n,Wen);
            rho_w2s=RGMM2st.par(1,1);
            beta_w2s=RGMM2st.par(2:end,1);
%             A2is=In/(In-rho_w2s*W_2);
            y2ps=rho_w2s*W_2*ystar+(x*beta_w2s);
            
            rho_w3s=RGMM3st.par(1,1);
            beta_w3s=RGMM3st.par(2:end,1);
%             A3is=In/(In-rho_w3s*W_3);
            y3ps=rho_w3s*W_3*ystar+(x*beta_w3s);
            
            X1=[x y2ps y3ps];
            Q1a=[h1 X1];

%     Augmented model estimation
    results1_b=gmm_sar_hetero_aug_p1_mc(ystar,X1,W_1,w1ys,Q1a,p1,para_ini,We);
%     Computation of the bootstraped statistic
    delt1_b=Con*[results1_b.par];
    varcov1_b=results1_b.var;
    Rv1Ri_b=eye(2)/(Con*varcov1_b*Con');
    stat1_b=delt1_b'*Rv1Ri_b*delt1_b;
    
%     Computation of the power of the 2 other statistics
            res_w2_star=reechantillonage(res_w2);

            ystar=A2i*(x*beta_w2+res_w2_star);
            w2ys=W_2*ystar;
            w1ys=W_1*ystar;
            w3ys=W_3*ystar;
            RGMM1st=gmm_sar_hetero_mc(ystar,x, W_1,w1ys,Q1,p1,par_ini_n,Wen);
            rho_w1s=RGMM1st.par(1,1);
            beta_w1s=RGMM1st.par(2:end,1);
            
%             A1is=In/(In-rho_w1s*W_1);
            y1ps=rho_w1s*W_1*ystar+(x*beta_w1s);
            
             RGMM3st=gmm_sar_hetero_mc(ystar,x, W_3,w3ys,Q3,p3,par_ini_n,Wen);
            rho_w3s=RGMM3st.par(1,1);
            beta_w3s=RGMM3st.par(2:end,1);
            
%             A3is=In/(In-rho_w3s*W_3);
            y3ps=rho_w3s*W_3*ystar+(x*beta_w3s);
            X2=[x y1ps y3ps];
            Q2a=[h2 X2];
              results2_b=gmm_sar_hetero_aug_p1_mc(ystar,X2,W_2,w2ys,Q2a,p2,para_ini,We);
                delt2_b=Con*[results2_b.par];
                varcov2_b=results2_b.var;
                Rv2Ri_b=eye(2)/(Con*varcov2_b*Con');
                stat2_b=delt2_b'*Rv2Ri_b*delt2_b;
                pow_2b(b,1)=stat2_b;
            %%%%%%
            res_w3_star=reechantillonage(res_w3);

            ystar=A3i*(x*beta_w3+res_w3_star);
            w3ys=W_3*ystar;
            w1ys=W_1*ystar;
            w2ys=W_2*ystar;
            RGMM1st=gmm_sar_hetero_mc(ystar,x, W_1,w1ys,Q1,p1,par_ini_n,Wen);
            rho_w1s=RGMM1st.par(1,1);
            beta_w1s=RGMM1st.par(2:end,1);
            
%             A1is=In/(In-rho_w1s*W_1);
            y1ps=rho_w1s*W_1*ystar+(x*beta_w1s);
            
             RGMM2st=gmm_sar_hetero_mc(ystar,x, W_2,w2ys,Q2,p2,par_ini_n,Wen);
            rho_w2s=RGMM2st.par(1,1);
            beta_w2s=RGMM2st.par(2:end,1);
            
%             A2is=In/(In-rho_w2s*W_2);
            y2ps=rho_w2s*W_2*ystar+(x*beta_w2s);
            X3=[x y1ps y2ps];
            Q3a=[h3 X3];
              results3_b=gmm_sar_hetero_aug_p1_mc(ystar,X3,W_3,w3ys,Q3a,p3,para_ini,We);
                delt3_b=Con*[results3_b.par];
                varcov3_b=results3_b.var;
                Rv3Ri_b=eye(2)/(Con*varcov3_b*Con');
                stat3_b=delt3_b'*Rv3Ri_b*delt3_b;
                pow_3b(b,1)=stat3_b;
               

   J1_boot(b,:)=[stat1_b, stat2,stat3];
end
min_b=min(J1_boot,[],2);
q_t=quantile(min_b,0.95);
q_s=quantile(J1_boot(:,1),0.95);

p_va=1-length(find(stat1>J1_boot(:,1)))/(B+1);
p_va_t=1-length(find(stat1>min_b))/(B+1);
p_va_pow_J2=1-length(find(stat2>pow_2b))/(B+1);
p_va_pow_J3=1-length(find(stat3>pow_3b))/(B+1);
signi.sca=[stat1 p_va p_va_t q_s q_t stat2 stat3 p_va_pow_J2 p_va_pow_J3];
signi.vec=[pow_2b pow_3b];
elseif ind==2 %Min is the stat with W_2
parfor b=1:B
    res_w2_star=reechantillonage(res_w2);
   
    ystar=A2i*(x*beta_w2+res_w2_star);
   w2ys=W_2*ystar;
   w1ys=W_1*ystar;
   w3ys=W_3*ystar;
            RGMM1st=gmm_sar_hetero_mc(ystar,x, W_1,w1ys,Q1,p1,par_ini_n,Wen);
            rho_w1s=RGMM1st.par(1,1);
            beta_w1s=RGMM1st.par(2:end,1);
            
%             A1is=In/(In-rho_w1s*W_1);
            y1ps=rho_w1s*W_1*ystar+(x*beta_w1s);
            
            RGMM3st=gmm_sar_hetero_mc(ystar,x, W_3,w3ys,Q3,p3,par_ini_n,Wen);
            rho_w3s=RGMM3st.par(1,1);
            beta_w3s=RGMM3st.par(2:end,1);
            
%             A3is=In/(In-rho_w3s*W_3);
            y3ps=rho_w3s*W_3*ystar+(x*beta_w3s);
            X2=[x y1ps y3ps];
            Q2a=[h2 X2];
              results2_b=gmm_sar_hetero_aug_p1_mc(ystar,X2,W_2,w2ys,Q2a,p2,para_ini,We);
                delt2_b=Con*[results2_b.par];
                varcov2_b=results2_b.var;
                Rv2Ri_b=eye(2)/(Con*varcov2_b*Con');
                stat2_b=delt2_b'*Rv2Ri_b*delt2_b;
                
% %                 Computation of the power of the 2 statistics
% FOR MODEL 3
            res_w3_star=reechantillonage(res_w3);
            ystar=A3i*(x*beta_w3+res_w3_star);
            w3ys=W_3*ystar;
            w1ys=W_1*ystar;
            w2ys=W_2*ystar;
            RGMM1st=gmm_sar_hetero_mc(ystar,x, W_1,w1ys,Q1,p1,par_ini_n,Wen);
            rho_w1s=RGMM1st.par(1,1);
            beta_w1s=RGMM1st.par(2:end,1);
            
%             A1is=In/(In-rho_w1s*W_1);
            y1ps=rho_w1s*W_1*ystar+(x*beta_w1s);
            
             RGMM2st=gmm_sar_hetero_mc(ystar,x, W_2,w2ys,Q2,p2,par_ini_n,Wen);
            rho_w2s=RGMM2st.par(1,1);
            beta_w2s=RGMM2st.par(2:end,1);
            
%             A2is=In/(In-rho_w2s*W_2);
            y2ps=rho_w2s*W_2*ystar+(x*beta_w2s);
            X3=[x y1ps y2ps];
            Q3a=[h3 X3];
              results3_b=gmm_sar_hetero_aug_p1_mc(ystar,X3,W_3,w3ys,Q3a,p3,para_ini,We);
                delt3_b=Con*[results3_b.par];
                varcov3_b=results3_b.var;
                Rv3Ri_b=eye(2)/(Con*varcov3_b*Con');
                stat3_b=delt3_b'*Rv3Ri_b*delt3_b;
                pow_3b(b,1)=stat3_b;
    
%          FOR MODEL 1       
                res_w1_star=reechantillonage(res_w1);
                ystar=A1i*(x*beta_w1+res_w1_star);
               w1ys=W_1*ystar;
               w2ys=W_2*ystar;
               w3ys=W_3*ystar;
               RGMM2st=gmm_sar_hetero_mc(ystar,x, W_2,w2ys,Q2,p2,par_ini_n,Wen);
               RGMM3st=gmm_sar_hetero_mc(ystar,x, W_3,w3ys,Q3,p3,par_ini_n,Wen);
                rho_w2s=RGMM2st.par(1,1);
            beta_w2s=RGMM2st.par(2:end,1);
%             A2is=In/(In-rho_w2s*W_2);
            y2ps=rho_w2s*W_2*ystar+(x*beta_w2s);
            
            rho_w3s=RGMM3st.par(1,1);
            beta_w3s=RGMM3st.par(2:end,1);
%             A3is=In/(In-rho_w3s*W_3);
            y3ps=rho_w3s*W_3*ystar+(x*beta_w3s);
            
            X1=[x y2ps y3ps];
            Q1a=[h1 X1];
    results1_b=gmm_sar_hetero_aug_p1_mc(ystar,X1,W_1,w1ys,Q1a,p1,para_ini,We);
%     Computation of the bootstraped statistic
    delt1_b=Con*[results1_b.par];
    varcov1_b=results1_b.var;
    Rv1Ri_b=eye(2)/(Con*varcov1_b*Con');
    stat1_b=delt1_b'*Rv1Ri_b*delt1_b;
    pow_1b(b,1)=stat1_b;
                
                
    J1_boot(b,:)=[stat1, stat2_b,stat3];
end
min_b=min(J1_boot,[],2);
q_t=quantile(min_b,0.95);
q_s=quantile(J1_boot(:,2),0.95);

p_va=1-length(find(stat2>J1_boot(:,2)))/(B+1);
p_va_t=1-length(find(stat2>min_b))/(B+1);
p_va_pow_J1=1-length(find(stat1>pow_1b))/(B+1);
p_va_pow_J3=1-length(find(stat3>pow_3b))/(B+1);
signi.sca=[stat2 p_va p_va_t q_s q_t stat1 stat3 p_va_pow_J1 p_va_pow_J3];
signi.vec=[pow_1b pow_3b];
% Case where the stat3 is the smallest of the 3 
else
    
parfor b=1:B
    res_w3_star=reechantillonage(res_w3);
   
    ystar=A3i*(x*beta_w3+res_w3_star);
     w3ys=W_3*ystar;
    w1ys=W_1*ystar;
    w2ys=W_2*ystar;
    RGMM2st=gmm_sar_hetero_mc(ystar,x, W_2,w2ys,Q2,p2,par_ini_n,Wen);
               RGMM1st=gmm_sar_hetero_mc(ystar,x, W_1,w1ys,Q1,p1,par_ini_n,Wen);
                rho_w2s=RGMM2st.par(1,1);
            beta_w2s=RGMM2st.par(2:end,1);
%             A2is=In/(In-rho_w2s*W_2);
            y2ps=rho_w2s*W_2*ystar+(x*beta_w2s);
            
            rho_w1s=RGMM1st.par(1,1);
            beta_w1s=RGMM1st.par(2:end,1);
%             A1is=In/(In-rho_w1s*W_1);
            y1ps=rho_w1s*W_1*ystar+(x*beta_w1s);
            
            X3=[x y2ps y1ps];
            Q3a=[h3 X3];

    
    results3_b=gmm_sar_hetero_aug_p1_mc(ystar,X3,W_3,w3ys,Q3a,p3,para_ini,We);
    
%     Computation of the bootstraped statistic
    delt3_b=Con*[results3_b.par];
    varcov3_b=results3_b.var;
    Rv3Ri_b=eye(2)/(Con*varcov3_b*Con');
    stat3_b=delt3_b'*Rv3Ri_b*delt3_b;
    
%     Computation of the power of the 2 other statistics
% For MODEL 1
 res_w1_star=reechantillonage(res_w1);
                ystar=A1i*(x*beta_w1+res_w1_star);
               w1ys=W_1*ystar;
               w2ys=W_2*ystar;
               w3ys=W_3*ystar;
               RGMM2st=gmm_sar_hetero_mc(ystar,x, W_2,w2ys,Q2,p2,par_ini_n,Wen);
               RGMM3st=gmm_sar_hetero_mc(ystar,x, W_3,w3ys,Q3,p3,par_ini_n,Wen);
                rho_w2s=RGMM2st.par(1,1);
            beta_w2s=RGMM2st.par(2:end,1);
%             A2is=In/(In-rho_w2s*W_2);
            y2ps=rho_w2s*W_2*ystar+(x*beta_w2s);
            
            rho_w3s=RGMM3st.par(1,1);
            beta_w3s=RGMM3st.par(2:end,1);
%             A3is=In/(In-rho_w3s*W_3);
            y3ps=rho_w3s*W_3*ystar+(x*beta_w3s);
            
            X1=[x y2ps y3ps];
            Q1a=[h1 X1];
    results1_b=gmm_sar_hetero_aug_p1_mc(ystar,X1,W_1,w1ys,Q1a,p1,para_ini,We);
%     Computation of the bootstraped statistic
    delt1_b=Con*[results1_b.par];
    varcov1_b=results1_b.var;
    Rv1Ri_b=eye(2)/(Con*varcov1_b*Con');
    stat1_b=delt1_b'*Rv1Ri_b*delt1_b;
    pow_1b(b,1)=stat1_b;
    
%     For Model2 
            res_w2_star=reechantillonage(res_w2);

            ystar=A2i*(x*beta_w2+res_w2_star);
            w2ys=W_2*ystar;
            w1ys=W_1*ystar;
            w3ys=W_3*ystar;
            RGMM1st=gmm_sar_hetero_mc(ystar,x, W_1,w1ys,Q1,p1,par_ini_n,Wen);
            rho_w1s=RGMM1st.par(1,1);
            beta_w1s=RGMM1st.par(2:end,1);
            
%             A1is=In/(In-rho_w1s*W_1);
            y1ps=rho_w1s*W_1*ystar+(x*beta_w1s);
            
             RGMM3st=gmm_sar_hetero_mc(ystar,x, W_3,w3ys,Q3,p3,par_ini_n,Wen);
            rho_w3s=RGMM3st.par(1,1);
            beta_w3s=RGMM3st.par(2:end,1);
            
%             A3is=In/(In-rho_w3s*W_3);
            y3ps=rho_w3s*W_3*ystar+(x*beta_w3s);
            X2=[x y1ps y3ps];
            Q2a=[h2 X2];
              results2_b=gmm_sar_hetero_aug_p1_mc(ystar,X2,W_2,w2ys,Q2a,p2,para_ini,We);
                delt2_b=Con*[results2_b.par];
                varcov2_b=results2_b.var;
                Rv2Ri_b=eye(2)/(Con*varcov2_b*Con');
                stat2_b=delt2_b'*Rv2Ri_b*delt2_b;
                pow_2b(b,1)=stat2_b;

                

    J1_boot(b,:)=[stat1, stat2, stat3_b];
end
min_b=min(J1_boot,[],2);
q_t=quantile(min_b,0.95);
q_s=quantile(J1_boot(:,3),0.95);

p_va=1-length(find(stat3>J1_boot(:,3)))/(B+1);
p_va_t=1-length(find(stat3>min_b))/(B+1);
p_va_pow_J1=1-length(find(stat1>pow_1b))/(B+1);
p_va_pow_J2=1-length(find(stat2>pow_2b))/(B+1);
signi.sca=[stat3 p_va p_va_t q_s q_t stat1 stat2 p_va_pow_J1 p_va_pow_J2]; 
signi.vec=[pow_1b pow_2b];
   
end

signi1.stat=signi.sca(:,1);
signi1.pval_J=signi.sca(:,2);
signi1.test=signi.sca(:,3);
signi1.rea=J1_boot; 
signi1.qu_s=signi.sca(:,4);
signi1.qu_t=signi.sca(:,5);
signi1.pval_M1alt=signi.sca(:,8);
signi1.pval_M2alt=signi.sca(:,9);
signi1.powerval_M1alt=signi.vec(:,1);
signi1.powerval_M2alt=signi.vec(:,2);

