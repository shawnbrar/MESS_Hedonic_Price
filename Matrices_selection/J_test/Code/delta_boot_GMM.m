function [signi_a,t_est] = delta_boot_GMM(y,x,W,Q,P1,crit)
%SUMMARY This fomction computes inference on impacts based on the Delta method.
% Also, I proposed a bootstrap refinement for small sample sizes.
% This code works for HETEROSKEDASTIC SAR models. Hence, the distribution
% used for the bootstrap is a Rademacher distribution, as advocated by
% Hageman (JoE 2012), who follows Davidson and Flachaire (JoE, 2008)

% Description
% y is the dependent variable
% x is the matrix of explanatory variables (including the constant)
% W is the interaction matrix
% Q is the set of linear moment conditions
% P is the quadratic moment
% crit is the critical value to compute the t-stat (1.65, 1.96, 2.58)
% Model estimation
est=gmmestimation_sar_he_old(y,x,W,Q,P1);
bhat = est.betah;
lhat=est.lambdah;
resh=est.resid;
covh=est.cov; % D'abord lmabda puis tous les beta)
% COmputation of impacts
[n,nvar]=size(x);
In=eye(n);
S_est=In/(In-lhat*W);
SWS=S_est*W*S_est; 
vS=vec(S_est);
C=[vec(SWS) vS]; 
if diag(covh)<=0
    warning('not all variances are positive');
end

vd=zeros(n^2,1);
for k = 1:nvar
cov_est=covh([1,k+1],[1,k+1]);
    for i =1:n^2
        vd(i,1)=[C(i,1)*bhat(k,1) C(i,2)]* cov_est*[C(i,1)*bhat(k,1) C(i,2)]';
    end
        sd_k=vd.^(0.5);
    t_est{k}=reshape((vS*bhat(k,1))./sd_k,n,n);
end

% The bootstrap needs to be thought further because there is a problem with
% it (nothing is significant). 
% % t_est represents the estimated t_stat
% % Let's now compute the bootstrapped statistic. 
% xb=x*bhat;
% vd_bo=zeros(n^2, 1);
% parfor b=1:B
%     res_star=reechantillonage(resh);
%     ystar=S_est*(xb+res_star);
%     boot=gmmestimation_sar_he(ystar,x, W,Q,P1);
%     lbo=boot.lambdah;
%     bbo=boot.betah;
%     cov_bo=boot.cov;
%     S_bo=In/(In-lbo*W);
%     SWS_bo=S_bo*W*S_bo;
%     vS_bo=vec(S_bo);
%     C_bo=[vec(SWS_bo), vS_bo]; 
%     for k = 1:nvar
%         vd_bo_k=zeros(n^2,1);
%         sd_bo_k=zeros(n^2,1);
%         cov_bou=cov_bo([1,k+1],[1,k+1]);
%         for i =1:n^2
%             vd_bo_k(i,1)=[C_bo(i,1)*bbo(k,1), C_bo(i,2)]* cov_bou*[C_bo(i,1)*bbo(k,1), C_bo(i,2)]';
%         end
%         sd_bo_k(:,b)=vd_bo_k.^(0.5);
%         t_bo(:,:,k,b)=reshape((vS_bo*bbo(k,1))./sd_bo_k(:,b),n,n);
%     end
% %     Add the estimated t-stat to te set of bootstrapped statistics 
% % For each explanatory variable k, we have an array of n x n x B+1 values,
% % i.e B+1 values for each t-
% end
%     for k = 1:nvar
%          t_bo(:,:,k,B+1)=t_est{k};
%     end
%     for k = 1:nvar
%         for j = 1:B+1
%              t_boo{k}(:,:,j)=t_bo(:,:,k,j);
%         end
%     end
%     
%   Computation of the acceptance zone for each impact, based on bootstrap.

  
  for k=1:nvar
    for i =1:n
      for j = 1:n
%           q_lo{k}(i,j)=quantile(t_boo{k}(i,j,:),0.025);
%           q_up{k}(i,j)=quantile(t_boo{k}(i,j,:),0.975);
%           signi_b{k}(i,j)=1-(q_lo{k}(i,j)<t_est{k}(i,j) && t_est{k}(i,j)<q_up{k}(i,j));
          signi_a{k}(i,j)=1-(-crit<t_est{k}(i,j) && t_est{k}(i,j)<crit);
      end
    end
  end

