function [estim, selectstat,n_param]=estimodel(y,x,z,x_z,W,initial_values,info,exo,endo_varx)

    fields = fieldnames(info);
    nf = length(fields);
    for i=1:nf
        if strcmp(fields{i},'cbiais')
            cbiais = info.cbiais;
        elseif strcmp(fields{i},'n')
            n = info.n;
        elseif strcmp(fields{i},'t')
            T = info.T;
        end  
    end
   
   % Estimation using Shi and Lee procedure 
   results = shi_lee_sspd_tvW(y,x,z,x_z,W,initial_values,info);
       
%     n=size(W,1);
%     t = size(y,1)/n-1; 
    nt=n*T;

    if cbiais==0
        c_par   = results.theta;
%         c_Sig   = results.SIGi;
%         c_Omg   = results.OMG;
        c_lik   = results.logL;
        c_tstat = results.tstat;
        c_std   = results.std;
%         c_Rsq_v = var(results.yhat)/var(results.yt);
%         c_Rsq_c = corr(results.yhat,results.yt)^2;
%         c_resid = results.resid;
%         c_residns = results.residns;
%         c_yhat = results.yhat;
    elseif cbiais==1
        c_par   = results.theta_bc;
%         c_Sig   = results.SIGi1;
%         c_Omg   = results.OMG1;
        c_lik   = results.logL_bc;
        c_tstat = results.tstat_bc;
        c_std   = results.std_bc;
%         c_Rsq_v = var(results.yhat1)/var(results.yt);
%         c_Rsq_c = corr(results.yhat1,results.yt)^2;
%         c_resid = results.resid1;
%         c_residns = results.residns1;
%         c_yhat = results.yhat1;
    end
    
%     gamma=c_par(1,1);
%     rho=c_par(2,1);
%     if sl==1
%         lambda=c_par(end-1,1);
%     end
%     beta = c_par(3:end-2,1);
    
    n_beta = exo;
    if isempty(endo_varx) == 1        
        n_param = [n_beta;n_beta;'WY';'Sigma_xi';'alpha';'delta']; 
    else
        n_param = [endo_varx;n_beta;'WY';'Sigma_xi';'alpha';'delta']; 
    end
          
    %sum up the results that will be printed
    estim = [c_par c_std c_tstat];
    aic   = -2*c_lik+2*size(c_par,1);
    bic   = -2*c_lik+size(c_par,1)*log(nt);
    
    selectstat=[c_lik;aic;bic];
%     goodfit = [c_Rsq_c;c_Rsq_v];
%     residuals = [c_resid, c_residns];
%     fitted = c_yhat;
    
 
    
    % Extract the common factors and factor loadings
     %res_fix_effect = felag_sdpd(y,x,W,c_par,info);
     
     % Alternative computation of yhat and residuals
%      yt = y(n+1:n+nt);
%      ytl = y(1:nt);
%      Wnt = kron(speye(t),W);
%      yhat_alt = gamma*ytl + rho*Wnt*ytl +lambda*Wnt*yt + x*beta + ...
%                                      repmat(res_fix_effect.fe,t,1);
%      resid_alt = yt - yhat_alt;
%      fitted_alt = yhat_alt;
%      
    
        
    
    
    
%   if strcmp(model,'no_unrt')==1 
%       %&& sum(strcmp(W_nam,{'defgdp','stability','socio_eco'}))==1
% 
 %       sptimeffect(n,x,W,c_par,V,W_nam,exo,info,n_param,effectpath);
% %         reseffect.Vnorm = (1/nt)*(c_Sig+c_Sig*c_Omg*c_Sig);,'defgdp'
% %     else 
% %         reseffect = 0; 
%   end
% %  reseffect = 0;
   

   
end
