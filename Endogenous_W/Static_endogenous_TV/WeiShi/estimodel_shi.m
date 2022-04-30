function [estim, selectstat,goodfit,nfactors,exitflag,results]=estimodel_shi(y,x,z,x_z,W,initial_values,info,exo,endo_varx,w_nam)%
% 
    fields = fieldnames(info);
    nf = length(fields);
    for i=1:nf
        if strcmp(fields{i},'cbiais')
            cbiais = info.cbiais;
        elseif strcmp(fields{i},'n')
            n = info.n;
        elseif strcmp(fields{i},'T')
            T = info.T;     
        end  
    end
   
    nt=n*T;
    
   % Estimation using Shi and Lee procedure 
%    Y =y;
%    Xy = x;
%    Xz = x_z;
%    Z = z;
%    TolX = 1e-9;

     results = shi_lee_sspd_tvW_correct(y,x,z,x_z,W,initial_values,info,w_nam);

       
%     n=size(W,1);
%     t = size(y,1)/n-1; 
    
    if cbiais == 0 
        c_par   = results.theta;
%         c_Sig   = results.SIGi;
%         c_Omg   = results.OMG;
        c_lik   = results.logL;
        c_liky   = results.lly;
        c_likz   = results.llz;
        c_tstat = results.tstat;
        c_std   = results.std;
        Rhat_y = results.Rhat_y;
        Rhat_z   = results.Rhat_z;
        exitflag = results.exitflag;
        c_Rsq_c = corr(reshape(results.yhat,size(results.yhat,1)*size(results.yhat,2),1),...
                        reshape(results.y,size(results.y,1)*size(results.y,2),1))^2;
%         c_Rsq_v = var(results.yhat)/var(results.yt);
%         c_Rsq_c = corr(results.yhat,results.yt)^2;
%         c_resid = results.resid;
%         c_residns = results.residns;
%         c_yhat = results.yhat;
    elseif cbiais == 1
        c_par   = results.theta_bc;
%         c_Sig   = results.SIGi1;
%         c_Omg   = results.OMG1;
        c_lik   = results.logL_bc;
        c_liky   = results.lly_bc;
        c_likz   = results.llz_bc;
        c_tstat = results.tstat_bc;
        c_std   = results.std_bc;
        
        Rhat_y = results.Rhat_y;
        Rhat_z   = results.Rhat_z;
        exitflag = results.exitflag;
       
        c_Rsq_c = corr(reshape(results.yhat_bc,size(results.yhat_bc,1)*size(results.yhat_bc,2),1),...
                        reshape(results.y,size(results.y,1)*size(results.y,2),1))^2;
        
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
    
%     n_beta = exo;
    
%     if isempty(endo_varx) == 1        
%         n_param = [n_beta;n_beta;'WY';'Sigma_xi';'alpha';'$\delta$']; 
%     else
%         n_param = [size(endo_varx,2);n_beta;'WY';'Sigma_xi';'alpha';'$\delta$']; 
%     end
          
    %sum up the results that will be printed
    estim = [c_par c_std c_tstat];
    aic   = -2*c_lik+2*size(c_par,1);
    bic   = -2*c_lik+size(c_par,1)*log(nt);
    
    selectstat=[c_lik;c_liky;c_likz;aic;bic];
    
    nfactors = [Rhat_y;Rhat_z];
    
    goodfit = [c_Rsq_c];

  
    
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
