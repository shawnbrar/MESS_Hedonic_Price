
function [output] = shi_lee_sspd_tvW_multZ(Y,X_Y,Z,X_Z,W,initial_values,info)

    fieldsinf = fieldnames(info);
    nf = length(fieldsinf);
    for i=1:nf
        if strcmp(fieldsinf{i},'n')
            n = info.n;
        elseif strcmp(fieldsinf{i},'T')
            T = info.T;
        elseif  strcmp(fieldsinf{i},'Ry')
            Ry = info.Ry;
        elseif strcmp(fieldsinf{i},'Rz')
            Rz = info.Rz; 
        elseif  strcmp(fieldsinf{i},'Ry')
            Ry = info.Ry;
        elseif strcmp(fieldsinf{i},'Rz')
            Rz = info.Rz; 
        elseif strcmp(fieldsinf{i},'factors_choice')
            factors_choice = info.factors_choice; 
        end
    end
    %X_Y = X;
    kz = (size(X_Z,2)/T)*size(X_Z,3);  %the number of explanatory variables in Z equation
    ky = size(X_Y,2)/T;                %the number of explanatory variables in y equation
    p  = size(Z,3);                    %the number of dependent variables in the Z equation
    J  = p + (p*(p-1))/2;              %the number of distinct elements in Sigma_epsilon

    
  % Optimization options          
   opt = optimoptions('fminunc','TolFun',1e-9,'TolX',1e-9,'MaxFunEvals',...
        100000,'MaxIter',100000,'Display','notify','Algorithm','quasi-newton');
    
   % Step 1: preliminary estimates, use a large number of factors as
   % specified in R_y and R_z
   
   
   for iv = 1:size(initial_values,2)
       theta = initial_values(:,iv);
       [t1,fval1,e1] = fminunc('Q_nT_multZ2',theta,opt,n,T,Rz,Ry,W,X_Y,Y,Z,X_Z,0);
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
       Rhat_y = determine_nfactors_multZ(est,W,Y,X_Y,Z,X_Z,'sar',info);

       % Endogenous W: z equation
       Rhat_z = determine_nfactors_multZ(est,W,Y,X_Y,Z,X_Z,'endo_W',info);
   elseif factors_choice == 0 % The user provides the right number of factors
       Rhat_z = info.Rz;
       Rhat_y = info.Ry;
   end
   
   % Step 3: Estimate the model with the estimated number of factors
   for iv = 1:size(initial_values,2)
       theta = initial_values(:,iv);
       [t1,fval1,e1] = fminunc('Q_nT_multZ2',theta,opt,n,T,Rhat_z,Rhat_y,W,X_Y,Y,Z,X_Z,0);
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
   
    [bc,sd,C,SIG] = bias_correction_multZ(est,n,T,Rhat_z,Rhat_y,W,X_Y,Y,Z,X_Z,0);
    [~,sd_bc,C_bc,SIG_bc] = bias_correction_multZ(bc,n,T,Rhat_z,Rhat_y,W,X_Y,Y,Z,X_Z,1);
 
    
    est(kz+ky+1) = 0.9*(exp(2*est(kz+ky+1))-1)/(exp(2*est(kz+ky+1))+1);
    est(kz+ky+2) = exp(est(kz+ky+2));
    est(kz+ky+3:kz+ky+2+J) = exp(est(kz+ky+3:kz+ky+2+J));
        
    % Compute the loglikelihood
%     logL = Q_nT(est,n,T,Rhat_z,Rhat_y,W,X_Y,Y,Z,X_Z,1);
%     logL_bc = Q_nT(bc,n,T,Rhat_z,Rhat_y,W,X_Y,Y,Z,X_Z,1);
    
    % Final: output
    output.n = n;
    output.T = T;
    output.Ry = Ry;
    output.Rz = Rz;
    output.Rhat_y = Rhat_y;
    output.Rhat_z = Rhat_z;
    output.theta = est;
    output.std = sd;
    %output.logL = logL;
    output.theta_bc = bc;
    output.std_bc = sd_bc;
   % output.logL_bc = logL_bc;
    output.C = C;
    output.SIG = SIG;
    output.C_bc = C_bc;
    output.SIG_bc = SIG_bc;
end