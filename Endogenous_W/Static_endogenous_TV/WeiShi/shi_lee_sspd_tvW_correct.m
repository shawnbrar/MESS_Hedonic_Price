
function [results] = shi_lee_sspd_tvW_correct(y,x,z,x_z,W,initial_values,info,w_nam)
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
        elseif strcmp(fieldsinf{i},'Rz')
            Rz = info.Rz; 
        elseif  strcmp(fieldsinf{i},'Ry')
            Ry = info.Ry;
        elseif strcmp(fieldsinf{i},'Rz')
            Rz = info.Rz; 
        elseif strcmp(fieldsinf{i},'TolX')
            TolX = info.TolX; 
        elseif strcmp(fieldsinf{i},'TolFun')
            TolFun = info.TolFun; 
        elseif strcmp(fieldsinf{i},'factors_choice')
            factors_choice = info.factors_choice;
        elseif strcmp(fieldsinf{i},'optimtyp')
            optimtyp = info.optimtyp; 
        elseif strcmp(fieldsinf{i},'boundtyp')
            boundtyp = info.boundtyp;  
        elseif strcmp(fieldsinf{i},'modeltyp')
            modeltyp = info.modeltyp;
        elseif strcmp(fieldsinf{i},'varset')
            varset = info.varset;
        end
    end
    
    kz = (size(x_z,2)/T)*size(x_z,3);  %the number of explanatory variables in Z equation
    ky = size(x,2)/T;                %the number of explanatory variables in y equation
    p  = size(z,3);                    %the number of dependent variables in the Z equation
    J  = p + (p*(p-1))/2;              %the number of distinct elements in Sigma_epsilon

    
  % Optimization options          
    if strcmp(optimtyp,'Shi_con') == 1 
        optcon = optimoptions('fmincon','Algorithm','sqp',...
            'TolFun',TolFun,'TolX',TolX,'MaxFunEvals',100000,'MaxIter',100000,...
            'Display','notify');
        if boundtyp == 1
            LB = [-Inf(kz,1);-Inf(ky,1);-0.99;0;0;-Inf(p,1)];
        elseif boundtyp == 2 
            LB = [-Inf(kz,1);-Inf(ky,1);-0.99;1e-7;1e-7;-Inf(p,1)];
        end
        UB = [Inf(kz,1);Inf(ky,1);0.99;Inf;Inf;Inf(p,1)]; 
    elseif strcmp(optimtyp,'Shi_unc') == 1
        optunc = optimset('TolFun',TolFun,'TolX',TolX,'MaxFunEvals',...
        100000,'MaxIter',100000,'Display','notify');
    end
   

   % Step 1: preliminary estimates, use a large number of factors as
   % specified in R_y and R_z
   %X_Y = X;
   for iv = 1:size(initial_values,2)
       theta = initial_values(:,iv);
       %theta = [zeros(kz,1);zeros(ky,1);-0.5;0.01;0.01;0];
   %    [t1,fval1,e1] = fminunc('Q_nT',theta,opt,n,T,Rz,Ry,W,x,y,z,x_z,1);
        if strcmp(optimtyp,'Shi_con') == 1 
            Q_nTcon_alt_oc=@(par) Q_nTcon_alt(par,n,T,Rz,Ry,W,x,y,z,x_z,info);
            [t1,fval1,e1] = fmincon(Q_nTcon_alt_oc,theta,[],[],[],[],LB,UB,[],...
                                       optcon);
        elseif strcmp(optimtyp,'Shi_unc') == 1 
        Q_obj_alt_oc=@(par) Q_obj_alt(par,y,W,x,x_z,z,Ry,Rz);
            [t1,fval1,e1] = fminunc(Q_obj_alt_oc,theta,optunc);
        end
        
%        if iv == 1 && e1>0
           est = t1;
           fval = fval1;
           exitflag = e1;
           ini_val = theta;
%        elseif iv~=1 && fval1<fval && e1>0
%            est = t1;
%            fval = fval1;
%            exitflag = e1;  
%            ini_val = theta;
%        end
   end
   
   % Step 2: Determine the number of factors
   if factors_choice == 1 % Go on to determine the number of factors
       % SAR equation
       Rhat_y = determine_nfactors(est,W,y,x,z,x_z,'sar',info);

       % Endogenous W: z equation
       Rhat_z = determine_nfactors(est,W,y,x,z,x_z,'endo_W',info);
   elseif factors_choice == 0 % The user provides the right number of factors
       Rhat_z = info.Rz;
       Rhat_y = info.Ry;
   end
   
   % Step 3: Estimate the model with the estimated number of factors
%    allt1 = [];
%    allfval1 = [];
%    alle1 = [];
%    foldnam = strcat('C:\Users\Sessi\Dropbox\Phd_thesis\shock_propagation_financial_inst\Output\Estimation\',...
%    modeltyp,'\Z_variables_as_ratio_to_tot_asset\inv_abs');  
%    cd(foldnam)
   for iv = 1:size(initial_values,2)
       theta = initial_values(:,iv);
        if strcmp(optimtyp,'Shi_con') == 1 
%             [t1_alt,fval1_alt,e1_alt,out1_alt] = fmincon('Q_nTcon',theta,[],[],[],[],LB,UB,[],...
%                                        optcon,n,T,Rhat_z,Rhat_y,W,x,y,z,x_z,info)
Q_nTcon_alt_oc=@(par) Q_nTcon_alt(par,n,T,Rhat_z,Rhat_y,W,x,y,z,x_z,info);
            [t1,fval1,e1] = fmincon(Q_nTcon_alt_oc,theta,[],[],[],[],LB,UB,[],...
                                       optcon);
        elseif strcmp(optimtyp,'Shi_unc') == 1 
            
%             [t1_alt,fval1_alt,e1_alt,out1_alt] = fminunc('Q_obj',theta,optunc,y,W,x,x_z,z,Rhat_y,Rhat_z)
Q_obj_alt_oc=@(par) Q_obj_alt(par,y,W,x,x_z,z,Rhat_y,Rhat_z);
            [t1,fval1,e1] = fminunc(Q_obj_alt_oc,theta,optunc);
        end
%        allt1 = [allt1,t1];
%        allfval1 = [allfval1,fval1];
%        alle1 = [alle1,e1];
       if iv == 1 && e1>0
           est = t1;
           fval = fval1;
           exitflag = e1;
           ini_val = theta;
       elseif iv~=1 && fval1<fval && e1>0
           est = t1;
           fval = fval1;
           exitflag = e1;  
           ini_val = theta;
       end
   end
%    filnam = strcat('sensitivity_',lower(strrep(optimtyp,'_','')),'_',...
%                      modeltyp,'_',varset,'.xlsx');
%    xlswrite(filnam,allt1,w_nam,strcat('B',num2str(2+size(theta,1)+4)))
%    xlswrite(filnam,allfval1,w_nam,strcat('B',num2str(2+2*size(theta,1)+7)))
%    xlswrite(filnam,alle1,w_nam,strcat('B',num2str(2+2*size(theta,1)+8))) 
   
   % Step 4: Bias correction procedure due to incidental paramters (factor
   % loadings and time factors)
    [bc,sd,gamma_z,f_z,gamma_y,f_y,C,SIG] = bias_correction(est,n,T,...
                                        Rhat_z,Rhat_y,W,x,y,z,x_z,info,0);
    [~,sd_bc,~,~,~,~,C_bc,SIG_bc] = bias_correction(bc,n,T,...
                                        Rhat_z,Rhat_y,W,x,y,z,x_z,info,1);

    %     est(kz+ky+1) = 0.9*(exp(2*est(kz+ky+1))-1)/(exp(2*est(kz+ky+1))+1);
%     est(kz+ky+2) = exp(est(kz+ky+2));
%     est(kz+ky+3:kz+ky+2+J) = exp(est(kz+ky+3:kz+ky+2+J));
        
    % Compute the loglikelihood
    [logL,lly,llz] = Q_nT_bis_alt(est,n,T,Rhat_z,Rhat_y,W,x,y,z,x_z,info,0);        
    [logL_bc,lly_bc,llz_bc] = Q_nT_bis_alt(bc,n,T,Rhat_z,Rhat_y,W,x,y,z,x_z,info,1);
    
      % Compute fitted values
%     yhat = (speye(nt) - p*Wnt)\(zt*bhat);
    beta_y = est(kz+1:ky+kz);
    lambda = est(ky+kz+1);
    lambda_std = sd(ky+kz+1);

    xsum = x*kron(beta_y,eye(T));
    G = zeros(n,T);
    for t = 1:T
        S = eye(n)-lambda*squeeze(W(:,:,t));
        G(:,t) = S\(xsum(:,t)+gamma_y*f_y(t,1:Rhat_y)');
    end  
    
    beta_y = bc(kz+1:ky+kz);
    lambda_bc = bc(ky+kz+1);
    lambda_bc_std = sd_bc(ky+kz+1);

    xsum = x*kron(beta_y,eye(T));
    Gbc = zeros(n,T);
    for t = 1:T
        S = eye(n)-lambda_bc*squeeze(W(:,:,t));
        Gbc(:,t) = S\(xsum(:,t)+gamma_y*f_y(t,1:Rhat_y)');
    end  

    % Final: output
    results.n = n;
    results.T = T;
    results.Ry = Ry;
    results.Rz = Rz;
    results.Rhat_y = Rhat_y;
    results.Rhat_z = Rhat_z;
    results.theta = est;
    results.std = sd;
    results.tstat = est./sd;
    results.logL = logL;
    results.lly = lly;
    results.llz = llz;
    results.theta_bc = bc;
    results.std_bc = sd_bc;
    results.tstat_bc = bc./sd_bc;
    results.logL_bc = logL_bc;
    results.lly_bc = lly_bc;
    results.llz_bc = llz_bc;
    results.C = C;
    results.SIG = SIG;
    results.C_bc = C_bc;
    results.SIG_bc = SIG_bc;
    results.Rhat_y = Rhat_y;
    results.gamma_y = gamma_y;
    results.f_y = f_y;
    results.Rhat_z = Rhat_z;
    results.gamma_z = gamma_z;
    results.f_z = f_z;
    % optim info
    results.exitflag = exitflag;
    results.ini_val = ini_val;
    results.yhat = G;
    results.yhat_bc = Gbc;
    results.y = y;
    results.lambda = lambda; 
    results.lambda_bc = lambda_bc; 
    results.lambda_std = lambda_std; 
    results.lambda_bc_std = lambda_bc_std; 

end