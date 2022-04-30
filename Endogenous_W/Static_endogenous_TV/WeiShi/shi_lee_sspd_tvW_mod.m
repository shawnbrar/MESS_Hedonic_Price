
function [results] = shi_lee_sspd_tvW_mod(y,x,z,x_z,W,initial_values,info)
% Updated on 17/10/17 to account for constrained optimization: line 30-36,
% 46-47,76-77 and 95 (opt)
    fieldsinf = fieldnames(info);
    nf = length(fieldsinf);
    for i=1:nf
        if strcmp(fieldsinf{i},'n')
            n = info.n;
        elseif strcmp(fieldsinf{i},'T')
            T = info.T;
        elseif  strcmp(fieldsinf{i},'Ry')
            Ry = info.Ry;
        elseif strcmp(fieldsinf{i},'factors_choice')
            factors_choice = info.factors_choice; 
        elseif strcmp(fieldsinf{i},'modeltyp')
            modeltyp = info.modeltyp;             
        end
    end
    
    ky = size(x,2)/T;                %the number of explanatory variables in y equation

    
  % Optimization options          
%    opt = optimoptions('fminunc','TolFun',1e-9,'TolX',1e-9,'MaxFunEvals',...
%         100000,'MaxIter',100000,'Display','notify','Algorithm','quasi-newton');
    optcon = optimoptions('fmincon','Algorithm','interior-point','TolFun',1e-9,'TolX',1e-9,'MaxFunEvals',...
         100000,'MaxIter',100000,'Display','notify');
     
     if strcmp(modeltyp,'exog') == 1
        LB = [-Inf(ky,1);-0.9;0];
        UB = [Inf(ky,1);0.9;Inf];         
     end

   % Step 1: preliminary estimates, use a large number of factors as
   % specified in R_y and R_z
   %X_Y = X;
   
   for iv = 1:size(initial_values,2)
       theta = initial_values(:,iv);
       %theta = [zeros(kz,1);zeros(ky,1);-0.5;0.01;0.01;0];
   %    [t1,fval1,e1] = fminunc('Q_nT',theta,opt,n,T,Rz,Ry,W,x,y,z,x_z,1);
       [t1,fval1,e1] = fmincon('Q_nTcon_exog',theta,[],[],[],[],LB,UB,[],optcon,n,T,Ry,W,x,y);
       if iv == 1 && e1>0
           est = t1;
           fval = fval1;
           exitflg = e1;
           ini_val = theta;
       elseif iv~=1 && fval1<fval && e1>0
           est = t1;
           fval = fval1;
           exitflag = e1;  
           ini_val = theta;
       end
   end
   
   % Step 2: Determine the number of factors
   if factors_choice == 1 % Go on to determine the number of factors
       % SAR equation
       Rhat_y = determine_nfactors_exog(est,W,y,x,info);
   elseif factors_choice == 0 % The user provides the right number of factors
       Rhat_y = info.Ry;
   end
   
   % Step 3: Estimate the model with the estimated number of factors
   for iv = 1:size(initial_values,2)
       theta = initial_values(:,iv);
       [t1,fval1,e1] = fmincon('Q_nTcon_exog',theta,[],[],[],[],LB,UB,[],optcon,n,T,Rhat_y,W,x,y);
       %[t1,fval1,e1] = fminunc('Q_nT',theta,opt,n,T,Rhat_z,Rhat_y,W,x,y,z,x_z,0);
       if iv == 1 && e1>0
           est = t1;
           fval = fval1;
           exitflg = e1;
           ini_val = theta;
       elseif iv~=1 && fval1<fval && e1>0
           est = t1;
           fval = fval1;
           exitflag = e1;  
           ini_val = theta;
       end
   end
   
   % Step 4: Bias correction procedure due to incidental paramters (factor
   % loadings and time factors)
   
    [bc,sd,C,SIG] = bias_correction_exog(est,n,T,Rhat_y,W,x,y);
    [~,sd_bc,C_bc,SIG_bc] = bias_correction_exog(bc,n,T,Rhat_y,W,x,y);
    
%     est(kz+ky+1) = 0.9*(exp(2*est(kz+ky+1))-1)/(exp(2*est(kz+ky+1))+1);
%     est(kz+ky+2) = exp(est(kz+ky+2));
%     est(kz+ky+3:kz+ky+2+J) = exp(est(kz+ky+3:kz+ky+2+J));
        
    % Compute the loglikelihood
    logL = Q_nTcon_exogbis(est,n,T,Rhat_y,W,x,y);
    logL_bc = Q_nTcon_exogbis(bc,n,T,Rhat_y,W,x,y);
    
    % Final: output
    results.n = n;
    results.T = T;
    results.Ry = Ry;
    results.Rhat_y = Rhat_y;
    results.theta = est;
    results.std = sd;
    results.tstat = est./sd;
    results.logL = logL;
    results.theta_bc = bc;
    results.std_bc = sd_bc;
    results.tstat_bc = bc./sd_bc;
    results.logL_bc = logL_bc;
    results.C = C;
    results.SIG = SIG;
    results.C_bc = C_bc;
    results.SIG_bc = SIG_bc;
    results.Rhat_y = Rhat_y;
end