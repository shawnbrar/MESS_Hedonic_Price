% Monte Carlo experiment for "A Spatial Panel Data Model with an Endogenous
%   Weights Matrix and Common Factors"
%
%   10/15/2015

tic

%% Set Monte Carlo design parameters
global n T theta_select Ry0 R Rz Ry W X Y Z SNR

Rz0 = 1;        % true number of factors 
rng(20151014);

%% Main script
name = strcat('MCn_',num2str(n),'T_',num2str(T),'R_',num2str(R),...
    'theta_',num2str(theta_select),'Ry0_',num2str(Ry0),'SNR_',num2str(SNR));
diary_name = strcat(name,'_logs.txt');
data_name = strcat(name,'_data.mat');
diary(diary_name)

% Defining different set of parameters
if theta_select == 1
    beta_z0  = 1;
    beta_y0  = 1;
    lambda0  = 0.2;
    sigma_v0 = sqrt(16/12*SNR);
    sigma_e0 = sqrt(16/12*SNR);
    rho0     = 0.2;
elseif theta_select == 2
    beta_z0  = 1;
    beta_y0  = 1;
    lambda0  = 0.6;
    sigma_v0 = sqrt(16/12*SNR);
    sigma_e0 = sqrt(16/12*SNR);
    rho0     = 0.2;
elseif theta_select == 3
    beta_z0  = 1;
    beta_y0  = 1;
    lambda0  = 0.2;
    sigma_v0 = sqrt(16/12*SNR);
    sigma_e0 = sqrt(16/12*SNR);
    rho0     = 0.6;
elseif theta_select == 4
    beta_z0  = 1;
    beta_y0  = 1;
    lambda0  = 0.6;
    sigma_v0 = sqrt(16/12*SNR);
    sigma_e0 = sqrt(16/12*SNR);
    rho0     = 0.6;    
end

%% Parameters
sigma_xi0 = sqrt((1-rho0^2)*sigma_v0^2);
alpha_xi0 = 1/sigma_xi0;   %  ????
alpha0 = 1/sigma_e0;       %  ????
delta0 = rho0*sigma_v0/sigma_e0;
m      = min(n,T);

%% Generate Data
NT      = n*T;
Gamma_y = rand(n,Ry0)*4-ones(n,Ry0)*2;
Gamma_z = rand(n,Rz0)*4-ones(n,Rz0)*2;
F_y     = rand(T,Ry0)*4-ones(T,Ry0)*2;
F_z     = F_y(:,1:Rz0);

X       = rand(n,T)*4-ones(n,T)*2 + (Gamma_y*F_y' + Gamma_y*ones(Ry0,1)*ones(1,T) + ones(n,1)*(F_y*ones(Ry0,1))')/3;
mu      = [0,0];
Sigma   = [sigma_v0^2,rho0*sigma_v0*sigma_e0;rho0*sigma_v0*sigma_e0,sigma_e0^2];

%% Create the spatial weights matrix
Wd  = zeros(n,n);
rr  = sqrt(n);
for i = 1:n
    for j = 1:n
        row_i       = ceil(i/rr);
        column_i    = i-(row_i-1)*rr;
        row_j       = ceil(j/rr);
        column_j    = j-(row_j-1)*rr;
        if abs(row_i-row_j)+abs(column_i-column_j) == 1
            Wd(i,j) = 1;
        end
    end
end

%% initial values
% Case of endogeneous weight matrix
theta1 = [0;0;0;0;0;0];
theta2 = [0;0;0.7;0;0;0];
theta3 = [0;0;0;0;0;0.7];
theta4 = [0;0;0.7;0;0;0.7];
theta5 = [0;0;0.2;0;0;0];
theta6 = [0;0;0;0;0;0.2];
theta7 = [0;0;0.2;0;0;0.2];
theta8 = randn(6,1);
theta9  = randn(6,1);
theta10 = randn(6,1);
theta11 = randn(6,1);
theta12 = randn(6,1);

% Case of exogeneous weight matrix
theta1_sar  = [0;0;0];
theta2_sar  = [0;0.7;0];
theta3_sar  = [0;0.2;0];
theta4_sar  = randn(3,1);
theta5_sar  = randn(3,1);
theta6_sar  = randn(3,1);
theta7_sar  = randn(3,1);
theta8_sar  = randn(3,1);
theta9_sar  = randn(3,1);
theta10_sar  = randn(3,1);
theta11_sar  = randn(3,1);
theta12_sar  = randn(3,1);

%% Monte Carlo Iterations
% saved results
factor_error_ER = zeros(R,2); % column 1: z equation, column 2: y equation
factor_error_GR = zeros(R,2);
theta           = zeros(R,6,2); % (.,.,1): QML estimates, (.,.,2): bias corrected estimates
theta_sd        = zeros(R,6,2);
theta_cp        = zeros(R,6,2);
theta_rf1       = zeros(R,6,2);
theta_rf1_sd    = zeros(R,6,2);
theta_rf1_cp    = zeros(R,6,2);
theta_rf2       = zeros(R,6,2);
theta_rf2_sd    = zeros(R,6,2);
theta_rf2_cp    = zeros(R,6,2);
theta_rf3       = zeros(R,6,2);
theta_rf3_sd    = zeros(R,6,2);
theta_rf3_cp    = zeros(R,6,2);
theta_ex        = zeros(R,3);

for mc = 1:R
    fprintf('Monte Carlo Iterations \t round %.0f \t %.2f%% \n',mc, (mc/R)*100);
    V   = zeros(n,T);
    E   = zeros(n,T);
    Y   = zeros(n,T);
    W   = zeros(n,n,T);
    
    Error = mvnrnd(mu,Sigma,NT);
    for t = 1:T
        V(:,t) = Error(1+(t-1)*n:t*n,1);
        E(:,t) = Error(1+(t-1)*n:t*n,2);
    end
    Z = X*beta_z0 + Gamma_z*F_z' + E;    
    
    for t = 1:T
        for i = 1:n
            for j = 1:n
                W(i,j,t) = Wd(i,j)*min(1/abs(Z(i,t)-Z(j,t)),2);
            end
        end
        W(:,:,t) = row_normalize(W(:,:,t));
    end
    
    for t = 1:T
        S0 = eye(n)-lambda0*W(:,:,t);
        Y(:,t) = S0\(X(:,t)*beta_y0 + Gamma_y*F_y(t,:)' + V(:,t));
    end
    
    opt = optimoptions('fminunc','TolFun',1e-9,'TolX',1e-9,'MaxFunEvals',...
        100000,'MaxIter',100000,'Display','notify','Algorithm','quasi-newton');
    
    % preliminary estimates, use a large number of factors
    Rz = 10;             % number of factors assumed in estimation
    Ry = 10;
    
    [t1,fval1,e1] = fminunc('Q_nT',theta1,opt);
    [t2,fval2,e2] = fminunc('Q_nT',theta2,opt);
    [t3,fval3,e3] = fminunc('Q_nT',theta3,opt);
    [t4,fval4,e4] = fminunc('Q_nT',theta4,opt);
    [t5,fval5,e5] = fminunc('Q_nT',theta5,opt);
    [t6,fval6,e6] = fminunc('Q_nT',theta6,opt);
    [t7,fval7,e7] = fminunc('Q_nT',theta7,opt);
    [t8,fval8,e8] = fminunc('Q_nT',theta8,opt);
    [t9,fval9,e9] = fminunc('Q_nT',theta9,opt);
    [t10,fval10,e10] = fminunc('Q_nT',theta10,opt);
    [t11,fval11,e11] = fminunc('Q_nT',theta11,opt);
    [t12,fval12,e12] = fminunc('Q_nT',theta12,opt);
    
    if fval1 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) ...
            && e1 > 0
        est = t1;
        val = 1;
        exitflag = e1;
    elseif fval2 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e2 > 0
        est = t2;
        val = 2;
        exitflag = e2;
    elseif fval3 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e3 > 0
        est = t3;
        val = 3;
        exitflag = e3;
    elseif fval4 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e4 > 0
        est = t4;
        val = 4;
        exitflag = e4;      
    elseif fval5 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e5 > 0
        est = t5;
        val = 5;
        exitflag = e5;    
    elseif fval6 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e6 > 0
        est = t6;
        val = 6;
        exitflag = e6;
    elseif fval7 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e7 > 0
        est = t7;
        val = 7;
        exitflag = e7;
    elseif fval8 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e8 > 0
        est = t8;
        val = 8;
        exitflag = e8;
    elseif fval9 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e9 > 0
        est = t9;
        val = 9;
        exitflag = e9;       
    elseif fval10 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e10 > 0
        est = t10;
        val = 10;
        exitflag = e10;    
    elseif fval11 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e11 > 0
        est = t11;
        val = 11;
        exitflag = e11;
    elseif fval12 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e12 > 0
        est = t12;
        val = 12;
        exitflag = e12;          
    end
    
    %% determine the number of factors using eigenvalue ratio
    fprintf('\t Factor number selection \t est = %.0f \t exitflag = %.0f \n',val,exitflag);
    beta_z     = est(1);
    beta_y     = est(2);
    lambda     = 0.9*(exp(2*est(3))-1)/(exp(2*est(3))+1);  % lambda < 0.9, > -0.9
    alpha_xi   = exp(est(4));   % sigma_xi^{-1}
    alpha      = exp(est(5));   % sigma_epsilon^{-1}
    delta      = est(6);

    sigma_eps_sq = 1/alpha^2;
    sigma_xi_sq  = 1/alpha_xi^2;

    D = (Z-X*beta_z)/sqrt(n*T*sigma_eps_sq);
    G = zeros(n,T);
    for t = 1:T
        S = eye(n)-lambda*W(:,:,t);
        G(:,t) = S*Y(:,t)-X(:,t)*beta_y-(Z(:,t)-X(:,t)*beta_z)*delta;
    end
    G = G/sqrt(n*T*sigma_xi_sq);

    % list of (sorted) eigenvalues
    [~,D1,~] = svd(D);
    D1 = diag(D1).*diag(D1);
    D1 = sort(D1);

    [~,G1,~] = svd(G);
    G1 = diag(G1).*diag(G1);
    G1 = sort(G1);

    %% Eigenvalue ratio criterion
    ER_z = zeros(21,1);
    ER_z(1) = sum(D1)/log(m)/D1(m);
    for k = 2:21
        ER_z(k) = D1(m+2-k)/D1(m+1-k);
    end
    [~,I]   = sort(ER_z);
    Rz_ER   = I(21)-1;

    ER_y = zeros(21,1);
    ER_y(1) = sum(G1)/log(m)/G1(m);
    for k = 2:21
        ER_y(k) = G1(m+2-k)/G1(m+1-k);
    end
    [~,I]   = sort(ER_y);
    Ry_ER   = I(21)-1;

    %% Growth ratio criterion
    GR_z = zeros(21,1);
    for k = 1:21
        if k == 1
            GR_z(k) = log((sum(D1)+sum(D1)/log(m))/sum(D1))/log(sum(D1)/sum(D1(1:m-1)));
        else
            GR_z(k) = log((sum(D1(1:m+2-k)))/sum(D1(1:m+1-k)))/log(sum(D1(1:m+1-k))/sum(D1(1:m-k)));
        end
    end
    [~,I]   = sort(GR_z);
    Rz_GR   = I(21)-1;

    GR_y = zeros(21,1);
    for k = 1:21
        if k == 1
            GR_y(k) = log((sum(G1)+sum(G1)/log(m))/sum(G1))/log(sum(G1)/sum(G1(1:m-1)));
        else
            GR_y(k) = log((sum(G1(1:m+2-k)))/sum(G1(1:m+1-k)))/log(sum(G1(1:m+1-k))/sum(G1(1:m-k)));
        end
    end
    [~,I]   = sort(GR_y);
    Ry_GR   = I(21)-1;
    
    if Rz_ER ~= Rz0
        factor_error_ER(mc,1) = 1;
    end
    if Ry_ER ~= Ry0
        factor_error_ER(mc,2) = 1;
    end
    if Rz_GR ~= Rz0
        factor_error_GR(mc,1) = 1;
    end
    if Ry_GR ~= Ry0
        factor_error_GR(mc,2) = 1;
    end
    
    %% Use the true number of factors
    fprintf('\t QML estimates using the true number of factors \n');
    Rz = Rz0;
    Ry = Ry0;
    
    [t1,fval1,e1] = fminunc('Q_nT',theta1,opt);
    [t2,fval2,e2] = fminunc('Q_nT',theta2,opt);
    [t3,fval3,e3] = fminunc('Q_nT',theta3,opt);
    [t4,fval4,e4] = fminunc('Q_nT',theta4,opt);
    [t5,fval5,e5] = fminunc('Q_nT',theta5,opt);
    [t6,fval6,e6] = fminunc('Q_nT',theta6,opt);
    [t7,fval7,e7] = fminunc('Q_nT',theta7,opt);
    [t8,fval8,e8] = fminunc('Q_nT',theta8,opt);
    [t9,fval9,e9] = fminunc('Q_nT',theta9,opt);
    [t10,fval10,e10] = fminunc('Q_nT',theta10,opt);
    [t11,fval11,e11] = fminunc('Q_nT',theta11,opt);
    [t12,fval12,e12] = fminunc('Q_nT',theta12,opt);
    
    if fval1 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e1 > 0
        est = t1;
        val = 1;
        exitflag = e1;
    elseif fval2 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e2 > 0
        est = t2;
        val = 2;
        exitflag = e2;
    elseif fval3 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e3 > 0
        est = t3;
        val = 3;
        exitflag = e3;
    elseif fval4 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e4 > 0
        est = t4;
        val = 4;
        exitflag = e4;      
    elseif fval5 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e5 > 0
        est = t5;
        val = 5;
        exitflag = e5;    
    elseif fval6 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e6 > 0
        est = t6;
        val = 6;
        exitflag = e6;
    elseif fval7 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e7 > 0
        est = t7;
        val = 7;
        exitflag = e7;
    elseif fval8 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e8 > 0
        est = t8;
        val = 8;
        exitflag = e8;
    elseif fval9 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e9 > 0
        est = t9;
        val = 9;
        exitflag = e9;       
    elseif fval10 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e10 > 0
        est = t10;
        val = 10;
        exitflag = e10;    
    elseif fval11 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e11 > 0
        est = t11;
        val = 11;
        exitflag = e11;
    elseif fval12 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e12 > 0
        est = t12;
        val = 12;
        exitflag = e12;          
    end
    
    [bc,sd] = bias_correction(est,0);
    [~,sd_bc] = bias_correction(bc,1);
    
    est(3) = 0.9*(exp(2*est(3))-1)/(exp(2*est(3))+1);
    est(4) = exp(est(4));
    est(5) = exp(est(5));
    
    theta(mc,:,1) = est';
    theta(mc,:,2) = bc';
    theta_sd(mc,:,1) = sd';
    theta_sd(mc,:,2) = sd_bc';
    
    % CP for the original QML estimator
    if (1-normcdf(abs((beta_z0-theta(mc,1,1))/theta_sd(mc,1,1))))*2 >= 0.05
        theta_cp(mc,1,1) = 1;
    end
    if (1-normcdf(abs((beta_y0-theta(mc,2,1))/theta_sd(mc,2,1))))*2 >= 0.05
        theta_cp(mc,2,1) = 1;
    end
    if (1-normcdf(abs((lambda0-theta(mc,3,1))/theta_sd(mc,3,1))))*2 >= 0.05
        theta_cp(mc,3,1) = 1;
    end
    if (1-normcdf(abs((alpha_xi0-theta(mc,4,1))/theta_sd(mc,4,1))))*2 >= 0.05
        theta_cp(mc,4,1) = 1;
    end
    if (1-normcdf(abs((alpha0-theta(mc,5,1))/theta_sd(mc,5,1))))*2 >= 0.05
        theta_cp(mc,5,1) = 1;
    end
    if (1-normcdf(abs((delta0-theta(mc,6,1))/theta_sd(mc,6,1))))*2 >= 0.05
        theta_cp(mc,6,1) = 1;
    end
 
    % CP for the bias corrected estimator
    if (1-normcdf(abs((beta_z0-theta(mc,1,2))/theta_sd(mc,1,2))))*2 >= 0.05
        theta_cp(mc,1,2) = 1;
    end
    if (1-normcdf(abs((beta_y0-theta(mc,2,2))/theta_sd(mc,2,2))))*2 >= 0.05
        theta_cp(mc,2,2) = 1;
    end
    if (1-normcdf(abs((lambda0-theta(mc,3,2))/theta_sd(mc,3,2))))*2 >= 0.05
        theta_cp(mc,3,2) = 1;
    end
    if (1-normcdf(abs((alpha_xi0-theta(mc,4,2))/theta_sd(mc,4,2))))*2 >= 0.05
        theta_cp(mc,4,2) = 1;
    end
    if (1-normcdf(abs((alpha0-theta(mc,5,2))/theta_sd(mc,5,2))))*2 >= 0.05
        theta_cp(mc,5,2) = 1;
    end
    if (1-normcdf(abs((delta0-theta(mc,6,2))/theta_sd(mc,6,2))))*2 >= 0.05
        theta_cp(mc,6,2) = 1;
    end       
    
    %% The number of factors in the y equation is over-specified by 1
    fprintf('\t The number of factors in the y equation is over-specified by 1... \n');
    Rz = Rz0;
    Ry = Ry0+1;
    
    [t1,fval1,e1] = fminunc('Q_nT',theta1,opt);
    [t2,fval2,e2] = fminunc('Q_nT',theta2,opt);
    [t3,fval3,e3] = fminunc('Q_nT',theta3,opt);
    [t4,fval4,e4] = fminunc('Q_nT',theta4,opt);
    [t5,fval5,e5] = fminunc('Q_nT',theta5,opt);
    [t6,fval6,e6] = fminunc('Q_nT',theta6,opt);
    [t7,fval7,e7] = fminunc('Q_nT',theta7,opt);
    [t8,fval8,e8] = fminunc('Q_nT',theta8,opt);
    [t9,fval9,e9] = fminunc('Q_nT',theta9,opt);
    [t10,fval10,e10] = fminunc('Q_nT',theta10,opt);
    [t11,fval11,e11] = fminunc('Q_nT',theta11,opt);
    [t12,fval12,e12] = fminunc('Q_nT',theta12,opt);
    
    if fval1 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e1 > 0
        est = t1;
        val = 1;
        exitflag = e1;
    elseif fval2 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e2 > 0
        est = t2;
        val = 2;
        exitflag = e2;
    elseif fval3 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e3 > 0
        est = t3;
        val = 3;
        exitflag = e3;
    elseif fval4 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e4 > 0
        est = t4;
        val = 4;
        exitflag = e4;      
    elseif fval5 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e5 > 0
        est = t5;
        val = 5;
        exitflag = e5;    
    elseif fval6 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e6 > 0
        est = t6;
        val = 6;
        exitflag = e6;
    elseif fval7 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e7 > 0
        est = t7;
        val = 7;
        exitflag = e7;
    elseif fval8 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e8 > 0
        est = t8;
        val = 8;
        exitflag = e8;
    elseif fval9 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e9 > 0
        est = t9;
        val = 9;
        exitflag = e9;       
    elseif fval10 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e10 > 0
        est = t10;
        val = 10;
        exitflag = e10;    
    elseif fval11 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e11 > 0
        est = t11;
        val = 11;
        exitflag = e11;
    elseif fval12 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e12 > 0
        est = t12;
        val = 12;
        exitflag = e12;          
    end
    
    [bc,sd] = bias_correction(est,0);
    [~,sd_bc] = bias_correction(bc,1);
    
    est(3) = 0.9*(exp(2*est(3))-1)/(exp(2*est(3))+1);
    est(4) = exp(est(4));
    est(5) = exp(est(5));
    
    theta_rf1(mc,:,1) = est';
    theta_rf1(mc,:,2) = bc';
    theta_rf1_sd(mc,:,1) = sd';
    theta_rf1_sd(mc,:,2) = sd_bc';    
    
    % CP for the original QML estimator
    if (1-normcdf(abs((beta_z0-theta_rf1(mc,1,1))/theta_rf1_sd(mc,1,1))))*2 >= 0.05
        theta_rf1_cp(mc,1,1) = 1;
    end
    if (1-normcdf(abs((beta_y0-theta_rf1(mc,2,1))/theta_rf1_sd(mc,2,1))))*2 >= 0.05
        theta_rf1_cp(mc,2,1) = 1;
    end
    if (1-normcdf(abs((lambda0-theta_rf1(mc,3,1))/theta_rf1_sd(mc,3,1))))*2 >= 0.05
        theta_rf1_cp(mc,3,1) = 1;
    end
    if (1-normcdf(abs((alpha_xi0-theta_rf1(mc,4,1))/theta_rf1_sd(mc,4,1))))*2 >= 0.05
        theta_rf1_cp(mc,4,1) = 1;
    end
    if (1-normcdf(abs((alpha0-theta_rf1(mc,5,1))/theta_rf1_sd(mc,5,1))))*2 >= 0.05
        theta_rf1_cp(mc,5,1) = 1;
    end
    if (1-normcdf(abs((delta0-theta_rf1(mc,6,1))/theta_rf1_sd(mc,6,1))))*2 >= 0.05
        theta_rf1_cp(mc,6,1) = 1;
    end
 
    % CP for the bias corrected estimator
    if (1-normcdf(abs((beta_z0-theta_rf1(mc,1,2))/theta_rf1_sd(mc,1,2))))*2 >= 0.05
        theta_rf1_cp(mc,1,2) = 1;
    end
    if (1-normcdf(abs((beta_y0-theta_rf1(mc,2,2))/theta_rf1_sd(mc,2,2))))*2 >= 0.05
        theta_rf1_cp(mc,2,2) = 1;
    end
    if (1-normcdf(abs((lambda0-theta_rf1(mc,3,2))/theta_rf1_sd(mc,3,2))))*2 >= 0.05
        theta_rf1_cp(mc,3,2) = 1;
    end
    if (1-normcdf(abs((alpha_xi0-theta_rf1(mc,4,2))/theta_rf1_sd(mc,4,2))))*2 >= 0.05
        theta_rf1_cp(mc,4,2) = 1;
    end
    if (1-normcdf(abs((alpha0-theta_rf1(mc,5,2))/theta_rf1_sd(mc,5,2))))*2 >= 0.05
        theta_rf1_cp(mc,5,2) = 1;
    end
    if (1-normcdf(abs((delta0-theta_rf1(mc,6,2))/theta_rf1_sd(mc,6,2))))*2 >= 0.05
        theta_rf1_cp(mc,6,2) = 1;
    end           
    
    %% The number of factors in the y equation is over-specified by 2
    fprintf('\t The number of factors in the y equation is over-specified by 2... \n');
    Rz = Rz0;
    Ry = Ry0+2;
    
    [t1,fval1,e1] = fminunc('Q_nT',theta1,opt);
    [t2,fval2,e2] = fminunc('Q_nT',theta2,opt);
    [t3,fval3,e3] = fminunc('Q_nT',theta3,opt);
    [t4,fval4,e4] = fminunc('Q_nT',theta4,opt);
    [t5,fval5,e5] = fminunc('Q_nT',theta5,opt);
    [t6,fval6,e6] = fminunc('Q_nT',theta6,opt);
    [t7,fval7,e7] = fminunc('Q_nT',theta7,opt);
    [t8,fval8,e8] = fminunc('Q_nT',theta8,opt);
    [t9,fval9,e9] = fminunc('Q_nT',theta9,opt);
    [t10,fval10,e10] = fminunc('Q_nT',theta10,opt);
    [t11,fval11,e11] = fminunc('Q_nT',theta11,opt);
    [t12,fval12,e12] = fminunc('Q_nT',theta12,opt);
    
    if fval1 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e1 > 0
        est = t1;
        val = 1;
        exitflag = e1;
    elseif fval2 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e2 > 0
        est = t2;
        val = 2;
        exitflag = e2;
    elseif fval3 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e3 > 0
        est = t3;
        val = 3;
        exitflag = e3;
    elseif fval4 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e4 > 0
        est = t4;
        val = 4;
        exitflag = e4;      
    elseif fval5 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e5 > 0
        est = t5;
        val = 5;
        exitflag = e5;    
    elseif fval6 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e6 > 0
        est = t6;
        val = 6;
        exitflag = e6;
    elseif fval7 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e7 > 0
        est = t7;
        val = 7;
        exitflag = e7;
    elseif fval8 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e8 > 0
        est = t8;
        val = 8;
        exitflag = e8;
    elseif fval9 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e9 > 0
        est = t9;
        val = 9;
        exitflag = e9;       
    elseif fval10 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e10 > 0
        est = t10;
        val = 10;
        exitflag = e10;    
    elseif fval11 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e11 > 0
        est = t11;
        val = 11;
        exitflag = e11;
    elseif fval12 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e12 > 0
        est = t12;
        val = 12;
        exitflag = e12;          
    end
    
    [bc,sd] = bias_correction(est,0);
    [~,sd_bc] = bias_correction(bc,1);
    
    est(3) = 0.9*(exp(2*est(3))-1)/(exp(2*est(3))+1);
    est(4) = exp(est(4));
    est(5) = exp(est(5));
    
    theta_rf2(mc,:,1) = est';
    theta_rf2(mc,:,2) = bc';
    theta_rf2_sd(mc,:,1) = sd';
    theta_rf2_sd(mc,:,2) = sd_bc';
    
    % CP for the original QML estimator
    if (1-normcdf(abs((beta_z0-theta_rf2(mc,1,1))/theta_rf2_sd(mc,1,1))))*2 >= 0.05
        theta_rf2_cp(mc,1,1) = 1;
    end
    if (1-normcdf(abs((beta_y0-theta_rf2(mc,2,1))/theta_rf2_sd(mc,2,1))))*2 >= 0.05
        theta_rf2_cp(mc,2,1) = 1;
    end
    if (1-normcdf(abs((lambda0-theta_rf2(mc,3,1))/theta_rf2_sd(mc,3,1))))*2 >= 0.05
        theta_rf2_cp(mc,3,1) = 1;
    end
    if (1-normcdf(abs((alpha_xi0-theta_rf2(mc,4,1))/theta_rf2_sd(mc,4,1))))*2 >= 0.05
        theta_rf2_cp(mc,4,1) = 1;
    end
    if (1-normcdf(abs((alpha0-theta_rf2(mc,5,1))/theta_rf2_sd(mc,5,1))))*2 >= 0.05
        theta_rf2_cp(mc,5,1) = 1;
    end
    if (1-normcdf(abs((delta0-theta_rf2(mc,6,1))/theta_rf2_sd(mc,6,1))))*2 >= 0.05
        theta_rf2_cp(mc,6,1) = 1;
    end
 
    % CP for the bias corrected estimator
    if (1-normcdf(abs((beta_z0-theta_rf2(mc,1,2))/theta_rf2_sd(mc,1,2))))*2 >= 0.05
        theta_rf2_cp(mc,1,2) = 1;
    end
    if (1-normcdf(abs((beta_y0-theta_rf2(mc,2,2))/theta_rf2_sd(mc,2,2))))*2 >= 0.05
        theta_rf2_cp(mc,2,2) = 1;
    end
    if (1-normcdf(abs((lambda0-theta_rf2(mc,3,2))/theta_rf2_sd(mc,3,2))))*2 >= 0.05
        theta_rf2_cp(mc,3,2) = 1;
    end
    if (1-normcdf(abs((alpha_xi0-theta_rf2(mc,4,2))/theta_rf2_sd(mc,4,2))))*2 >= 0.05
        theta_rf2_cp(mc,4,2) = 1;
    end
    if (1-normcdf(abs((alpha0-theta_rf2(mc,5,2))/theta_rf2_sd(mc,5,2))))*2 >= 0.05
        theta_rf2_cp(mc,5,2) = 1;
    end
    if (1-normcdf(abs((delta0-theta_rf2(mc,6,2))/theta_rf2_sd(mc,6,2))))*2 >= 0.05
        theta_rf2_cp(mc,6,2) = 1;
    end  
    
    %% The number of factors in the y equation is over-specified by 3
    fprintf('\t The number of factors in the y equation is over-specified by 3... \n');
    Rz = Rz0;
    Ry = Ry0+3;
    
    [t1,fval1,e1] = fminunc('Q_nT',theta1,opt);
    [t2,fval2,e2] = fminunc('Q_nT',theta2,opt);
    [t3,fval3,e3] = fminunc('Q_nT',theta3,opt);
    [t4,fval4,e4] = fminunc('Q_nT',theta4,opt);
    [t5,fval5,e5] = fminunc('Q_nT',theta5,opt);
    [t6,fval6,e6] = fminunc('Q_nT',theta6,opt);
    [t7,fval7,e7] = fminunc('Q_nT',theta7,opt);
    [t8,fval8,e8] = fminunc('Q_nT',theta8,opt);
    [t9,fval9,e9] = fminunc('Q_nT',theta9,opt);
    [t10,fval10,e10] = fminunc('Q_nT',theta10,opt);
    [t11,fval11,e11] = fminunc('Q_nT',theta11,opt);
    [t12,fval12,e12] = fminunc('Q_nT',theta12,opt);
    
    if fval1 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e1 > 0
        est = t1;
        val = 1;
        exitflag = e1;
    elseif fval2 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e2 > 0
        est = t2;
        val = 2;
        exitflag = e2;
    elseif fval3 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e3 > 0
        est = t3;
        val = 3;
        exitflag = e3;
    elseif fval4 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e4 > 0
        est = t4;
        val = 4;
        exitflag = e4;      
    elseif fval5 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e5 > 0
        est = t5;
        val = 5;
        exitflag = e5;    
    elseif fval6 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e6 > 0
        est = t6;
        val = 6;
        exitflag = e6;
    elseif fval7 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e7 > 0
        est = t7;
        val = 7;
        exitflag = e7;
    elseif fval8 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e8 > 0
        est = t8;
        val = 8;
        exitflag = e8;
    elseif fval9 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e9 > 0
        est = t9;
        val = 9;
        exitflag = e9;       
    elseif fval10 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e10 > 0
        est = t10;
        val = 10;
        exitflag = e10;    
    elseif fval11 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e11 > 0
        est = t11;
        val = 11;
        exitflag = e11;
    elseif fval12 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e12 > 0
        est = t12;
        val = 12;
        exitflag = e12;          
    end
    
    [bc,sd] = bias_correction(est,0);
    [~,sd_bc] = bias_correction(bc,1);
    
    est(3) = 0.9*(exp(2*est(3))-1)/(exp(2*est(3))+1);
    est(4) = exp(est(4));
    est(5) = exp(est(5));
    
    theta_rf3(mc,:,1) = est';
    theta_rf3(mc,:,2) = bc';
    theta_rf3_sd(mc,:,1) = sd';
    theta_rf3_sd(mc,:,2) = sd_bc';
    
    % CP for the original QML estimator
    if (1-normcdf(abs((beta_z0-theta_rf3(mc,1,1))/theta_rf3_sd(mc,1,1))))*2 >= 0.05
        theta_rf3_cp(mc,1,1) = 1;
    end
    if (1-normcdf(abs((beta_y0-theta_rf3(mc,2,1))/theta_rf3_sd(mc,2,1))))*2 >= 0.05
        theta_rf3_cp(mc,2,1) = 1;
    end
    if (1-normcdf(abs((lambda0-theta_rf3(mc,3,1))/theta_rf3_sd(mc,3,1))))*2 >= 0.05
        theta_rf3_cp(mc,3,1) = 1;
    end
    if (1-normcdf(abs((alpha_xi0-theta_rf3(mc,4,1))/theta_rf3_sd(mc,4,1))))*2 >= 0.05
        theta_rf3_cp(mc,4,1) = 1;
    end
    if (1-normcdf(abs((alpha0-theta_rf3(mc,5,1))/theta_rf3_sd(mc,5,1))))*2 >= 0.05
        theta_rf3_cp(mc,5,1) = 1;
    end
    if (1-normcdf(abs((delta0-theta_rf3(mc,6,1))/theta_rf3_sd(mc,6,1))))*2 >= 0.05
        theta_rf3_cp(mc,6,1) = 1;
    end
 
    % CP for the bias corrected estimator
    if (1-normcdf(abs((beta_z0-theta_rf3(mc,1,2))/theta_rf3_sd(mc,1,2))))*2 >= 0.05
        theta_rf3_cp(mc,1,2) = 1;
    end
    if (1-normcdf(abs((beta_y0-theta_rf3(mc,2,2))/theta_rf3_sd(mc,2,2))))*2 >= 0.05
        theta_rf3_cp(mc,2,2) = 1;
    end
    if (1-normcdf(abs((lambda0-theta_rf3(mc,3,2))/theta_rf3_sd(mc,3,2))))*2 >= 0.05
        theta_rf3_cp(mc,3,2) = 1;
    end
    if (1-normcdf(abs((alpha_xi0-theta_rf3(mc,4,2))/theta_rf3_sd(mc,4,2))))*2 >= 0.05
        theta_rf3_cp(mc,4,2) = 1;
    end
    if (1-normcdf(abs((alpha0-theta_rf3(mc,5,2))/theta_rf3_sd(mc,5,2))))*2 >= 0.05
        theta_rf3_cp(mc,5,2) = 1;
    end
    if (1-normcdf(abs((delta0-theta_rf3(mc,6,2))/theta_rf3_sd(mc,6,2))))*2 >= 0.05
        theta_rf3_cp(mc,6,2) = 1;
    end
    
    %% Estimation under misspecification: assuming that the spatial weights matrices is exogenous, factor is still considered
    fprintf('\t Assuming exogenous spatial weights matrices ... \n');
    Ry = Ry0;
    
    [t1,fval1,e1] = fminunc('Q_nT_SAR',theta1_sar,opt);
    [t2,fval2,e2] = fminunc('Q_nT_SAR',theta2_sar,opt);
    [t3,fval3,e3] = fminunc('Q_nT_SAR',theta3_sar,opt);
    [t4,fval4,e4] = fminunc('Q_nT_SAR',theta4_sar,opt);
    [t5,fval5,e5] = fminunc('Q_nT_SAR',theta5_sar,opt);
    [t6,fval6,e6] = fminunc('Q_nT_SAR',theta6_sar,opt);
    [t7,fval7,e7] = fminunc('Q_nT_SAR',theta7_sar,opt);
    [t8,fval8,e8] = fminunc('Q_nT_SAR',theta8_sar,opt);
    [t9,fval9,e9] = fminunc('Q_nT_SAR',theta9_sar,opt);
    [t10,fval10,e10] = fminunc('Q_nT_SAR',theta10_sar,opt);
    [t11,fval11,e11] = fminunc('Q_nT_SAR',theta11_sar,opt);
    [t12,fval12,e12] = fminunc('Q_nT_SAR',theta12_sar,opt);
    
    if fval1 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e1 > 0
        est = t1;
        val = 1;
        exitflag = e1;
    elseif fval2 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e2 > 0
        est = t2;
        val = 2;
        exitflag = e2;
    elseif fval3 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e3 > 0
        est = t3;
        val = 3;
        exitflag = e3;
    elseif fval4 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e4 > 0
        est = t4;
        val = 4;
        exitflag = e4;      
    elseif fval5 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e5 > 0
        est = t5;
        val = 5;
        exitflag = e5;    
    elseif fval6 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e6 > 0
        est = t6;
        val = 6;
        exitflag = e6;
    elseif fval7 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e7 > 0
        est = t7;
        val = 7;
        exitflag = e7;
    elseif fval8 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e8 > 0
        est = t8;
        val = 8;
        exitflag = e8;
    elseif fval9 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e9 > 0
        est = t9;
        val = 9;
        exitflag = e9;       
    elseif fval10 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e10 > 0
        est = t10;
        val = 10;
        exitflag = e10;    
    elseif fval11 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e11 > 0
        est = t11;
        val = 11;
        exitflag = e11;
    elseif fval12 == min([fval1 fval2 fval3 fval4 fval5 fval6 fval7 fval8 fval9 fval10 fval11 fval12]) && e12 > 0
        est = t12;
        val = 12;
        exitflag = e12;          
    end
    
    est(2) = 0.9*(exp(2*est(2))-1)/(exp(2*est(2))+1);   % lambda
    est(3) = exp(est(3));                               % alpha_xi
    
    theta_ex(mc,:) = est'; 
end

fprintf('\n \t n = %.0f T = %.0f \t Monte Carlo R = %.0f \n',n,T,R);
fprintf('\n Table 1 \t QML estimates \n');
fprintf('\t beta_z      \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',beta_z0,    mean(theta(:,1,1))-beta_z0,      std(theta(:,1,1)), sqrt(((theta(:,1,1)-ones(R,1)*beta_z0)'*(theta(:,1,1)-ones(R,1)*beta_z0))/R),      mean(theta_sd(:,1,1)), mean(theta_cp(:,1,1)));
fprintf('\t beta_y      \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',beta_y0,    mean(theta(:,2,1))-beta_y0,      std(theta(:,2,1)), sqrt(((theta(:,2,1)-ones(R,1)*beta_y0)'*(theta(:,2,1)-ones(R,1)*beta_y0))/R),      mean(theta_sd(:,2,1)), mean(theta_cp(:,2,1)));
fprintf('\t lambda      \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',lambda0,    mean(theta(:,3,1))-lambda0,      std(theta(:,3,1)), sqrt(((theta(:,3,1)-ones(R,1)*lambda0)'*(theta(:,3,1)-ones(R,1)*lambda0))/R),      mean(theta_sd(:,3,1)), mean(theta_cp(:,3,1)));
fprintf('\t alpha_xi    \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',alpha_xi0,  mean(theta(:,4,1))-alpha_xi0,    std(theta(:,4,1)), sqrt(((theta(:,4,1)-ones(R,1)*alpha_xi0)'*(theta(:,4,1)-ones(R,1)*alpha_xi0))/R),  mean(theta_sd(:,4,1)), mean(theta_cp(:,4,1)));
fprintf('\t alpha       \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',alpha0,     mean(theta(:,5,1))-alpha0,       std(theta(:,5,1)), sqrt(((theta(:,5,1)-ones(R,1)*alpha0)'*(theta(:,5,1)-ones(R,1)*alpha0))/R),        mean(theta_sd(:,5,1)), mean(theta_cp(:,5,1)));
fprintf('\t delta       \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',delta0,     mean(theta(:,6,1))-delta0,       std(theta(:,6,1)), sqrt(((theta(:,6,1)-ones(R,1)*delta0)'*(theta(:,6,1)-ones(R,1)*delta0))/R),        mean(theta_sd(:,6,1)), mean(theta_cp(:,6,1)));

fprintf('\n Table 2 \t Bias corrected estimates \n');
fprintf('\t beta_z      \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',beta_z0,    mean(theta(:,1,2))-beta_z0,      std(theta(:,1,2)), sqrt(((theta(:,1,2)-ones(R,1)*beta_z0)'*(theta(:,1,2)-ones(R,1)*beta_z0))/R),      mean(theta_sd(:,1,2)), mean(theta_cp(:,1,2)));
fprintf('\t beta_y      \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',beta_y0,    mean(theta(:,2,2))-beta_y0,      std(theta(:,2,2)), sqrt(((theta(:,2,2)-ones(R,1)*beta_y0)'*(theta(:,2,2)-ones(R,1)*beta_y0))/R),      mean(theta_sd(:,2,2)), mean(theta_cp(:,2,2)));
fprintf('\t lambda      \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',lambda0,    mean(theta(:,3,2))-lambda0,      std(theta(:,3,2)), sqrt(((theta(:,3,2)-ones(R,1)*lambda0)'*(theta(:,3,2)-ones(R,1)*lambda0))/R),      mean(theta_sd(:,3,2)), mean(theta_cp(:,3,2)));
fprintf('\t alpha_xi    \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',alpha_xi0,  mean(theta(:,4,2))-alpha_xi0,    std(theta(:,4,2)), sqrt(((theta(:,4,2)-ones(R,1)*alpha_xi0)'*(theta(:,4,2)-ones(R,1)*alpha_xi0))/R),  mean(theta_sd(:,4,2)), mean(theta_cp(:,4,2)));
fprintf('\t alpha       \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',alpha0,     mean(theta(:,5,2))-alpha0,       std(theta(:,5,2)), sqrt(((theta(:,5,2)-ones(R,1)*alpha0)'*(theta(:,5,2)-ones(R,1)*alpha0))/R),        mean(theta_sd(:,5,2)), mean(theta_cp(:,5,2)));
fprintf('\t delta       \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',delta0,     mean(theta(:,6,2))-delta0,       std(theta(:,6,2)), sqrt(((theta(:,6,2)-ones(R,1)*delta0)'*(theta(:,6,2)-ones(R,1)*delta0))/R),        mean(theta_sd(:,6,2)), mean(theta_cp(:,6,2)));

fprintf('\n Table 3 \t Factor number selection accuracy \n');
fprintf('\t Rz0 = %.0f, \t Ry0 = %.0f \n',Rz0,Ry0);
fprintf('\t Rz error frequency \t ER = %.5f \t GR = %.5f \n',mean(factor_error_ER(:,1)),mean(factor_error_GR(:,1)));
fprintf('\t Ry error frequency \t ER = %.5f \t GR = %.5f \n',mean(factor_error_ER(:,2)),mean(factor_error_GR(:,2)));

fprintf('\n Table 4A \t QML estimates, number of factors in the y equation is over-specified by 1 \n');
fprintf('\t beta_z      \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',beta_z0,    mean(theta_rf1(:,1,1))-beta_z0,      std(theta_rf1(:,1,1)), sqrt(((theta_rf1(:,1,1)-ones(R,1)*beta_z0)'*(theta_rf1(:,1,1)-ones(R,1)*beta_z0))/R),      mean(theta_rf1_sd(:,1,1)), mean(theta_rf1_cp(:,1,1)));
fprintf('\t beta_y      \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',beta_y0,    mean(theta_rf1(:,2,1))-beta_y0,      std(theta_rf1(:,2,1)), sqrt(((theta_rf1(:,2,1)-ones(R,1)*beta_y0)'*(theta_rf1(:,2,1)-ones(R,1)*beta_y0))/R),      mean(theta_rf1_sd(:,2,1)), mean(theta_rf1_cp(:,2,1)));
fprintf('\t lambda      \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',lambda0,    mean(theta_rf1(:,3,1))-lambda0,      std(theta_rf1(:,3,1)), sqrt(((theta_rf1(:,3,1)-ones(R,1)*lambda0)'*(theta_rf1(:,3,1)-ones(R,1)*lambda0))/R),      mean(theta_rf1_sd(:,3,1)), mean(theta_rf1_cp(:,3,1)));
fprintf('\t alpha_xi    \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',alpha_xi0,  mean(theta_rf1(:,4,1))-alpha_xi0,    std(theta_rf1(:,4,1)), sqrt(((theta_rf1(:,4,1)-ones(R,1)*alpha_xi0)'*(theta_rf1(:,4,1)-ones(R,1)*alpha_xi0))/R),  mean(theta_rf1_sd(:,4,1)), mean(theta_rf1_cp(:,4,1)));
fprintf('\t alpha       \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',alpha0,     mean(theta_rf1(:,5,1))-alpha0,       std(theta_rf1(:,5,1)), sqrt(((theta_rf1(:,5,1)-ones(R,1)*alpha0)'*(theta_rf1(:,5,1)-ones(R,1)*alpha0))/R),        mean(theta_rf1_sd(:,5,1)), mean(theta_rf1_cp(:,5,1)));
fprintf('\t delta       \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',delta0,     mean(theta_rf1(:,6,1))-delta0,       std(theta_rf1(:,6,1)), sqrt(((theta_rf1(:,6,1)-ones(R,1)*delta0)'*(theta_rf1(:,6,1)-ones(R,1)*delta0))/R),        mean(theta_rf1_sd(:,6,1)), mean(theta_rf1_cp(:,6,1)));

fprintf('\n Table 4B \t Bias corrected estimates, number of factors in the y equation is over-specified by 1 \n');
fprintf('\t beta_z      \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',beta_z0,    mean(theta_rf1(:,1,2))-beta_z0,      std(theta_rf1(:,1,2)), sqrt(((theta_rf1(:,1,2)-ones(R,1)*beta_z0)'*(theta_rf1(:,1,2)-ones(R,1)*beta_z0))/R),      mean(theta_rf1_sd(:,1,2)), mean(theta_rf1_cp(:,1,2)));
fprintf('\t beta_y      \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',beta_y0,    mean(theta_rf1(:,2,2))-beta_y0,      std(theta_rf1(:,2,2)), sqrt(((theta_rf1(:,2,2)-ones(R,1)*beta_y0)'*(theta_rf1(:,2,2)-ones(R,1)*beta_y0))/R),      mean(theta_rf1_sd(:,2,2)), mean(theta_rf1_cp(:,2,2)));
fprintf('\t lambda      \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',lambda0,    mean(theta_rf1(:,3,2))-lambda0,      std(theta_rf1(:,3,2)), sqrt(((theta_rf1(:,3,2)-ones(R,1)*lambda0)'*(theta_rf1(:,3,2)-ones(R,1)*lambda0))/R),      mean(theta_rf1_sd(:,3,2)), mean(theta_rf1_cp(:,3,2)));
fprintf('\t alpha_xi    \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',alpha_xi0,  mean(theta_rf1(:,4,2))-alpha_xi0,    std(theta_rf1(:,4,2)), sqrt(((theta_rf1(:,4,2)-ones(R,1)*alpha_xi0)'*(theta_rf1(:,4,2)-ones(R,1)*alpha_xi0))/R),  mean(theta_rf1_sd(:,4,2)), mean(theta_rf1_cp(:,4,2)));
fprintf('\t alpha       \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',alpha0,     mean(theta_rf1(:,5,2))-alpha0,       std(theta_rf1(:,5,2)), sqrt(((theta_rf1(:,5,2)-ones(R,1)*alpha0)'*(theta_rf1(:,5,2)-ones(R,1)*alpha0))/R),        mean(theta_rf1_sd(:,5,2)), mean(theta_rf1_cp(:,5,2)));
fprintf('\t delta       \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',delta0,     mean(theta_rf1(:,6,2))-delta0,       std(theta_rf1(:,6,2)), sqrt(((theta_rf1(:,6,2)-ones(R,1)*delta0)'*(theta_rf1(:,6,2)-ones(R,1)*delta0))/R),        mean(theta_rf1_sd(:,6,2)), mean(theta_rf1_cp(:,6,2)));

fprintf('\n Table 5A \t QML estimates, number of factors in the y equation is over-specified by 2 \n');
fprintf('\t beta_z      \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',beta_z0,    mean(theta_rf2(:,1,1))-beta_z0,      std(theta_rf2(:,1,1)), sqrt(((theta_rf2(:,1,1)-ones(R,1)*beta_z0)'*(theta_rf2(:,1,1)-ones(R,1)*beta_z0))/R),      mean(theta_rf2_sd(:,1,1)), mean(theta_rf2_cp(:,1,1)));
fprintf('\t beta_y      \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',beta_y0,    mean(theta_rf2(:,2,1))-beta_y0,      std(theta_rf2(:,2,1)), sqrt(((theta_rf2(:,2,1)-ones(R,1)*beta_y0)'*(theta_rf2(:,2,1)-ones(R,1)*beta_y0))/R),      mean(theta_rf2_sd(:,2,1)), mean(theta_rf2_cp(:,2,1)));
fprintf('\t lambda      \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',lambda0,    mean(theta_rf2(:,3,1))-lambda0,      std(theta_rf2(:,3,1)), sqrt(((theta_rf2(:,3,1)-ones(R,1)*lambda0)'*(theta_rf2(:,3,1)-ones(R,1)*lambda0))/R),      mean(theta_rf2_sd(:,3,1)), mean(theta_rf2_cp(:,3,1)));
fprintf('\t alpha_xi    \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',alpha_xi0,  mean(theta_rf2(:,4,1))-alpha_xi0,    std(theta_rf2(:,4,1)), sqrt(((theta_rf2(:,4,1)-ones(R,1)*alpha_xi0)'*(theta_rf2(:,4,1)-ones(R,1)*alpha_xi0))/R),  mean(theta_rf2_sd(:,4,1)), mean(theta_rf2_cp(:,4,1)));
fprintf('\t alpha       \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',alpha0,     mean(theta_rf2(:,5,1))-alpha0,       std(theta_rf2(:,5,1)), sqrt(((theta_rf2(:,5,1)-ones(R,1)*alpha0)'*(theta_rf2(:,5,1)-ones(R,1)*alpha0))/R),        mean(theta_rf2_sd(:,5,1)), mean(theta_rf2_cp(:,5,1)));
fprintf('\t delta       \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',delta0,     mean(theta_rf2(:,6,1))-delta0,       std(theta_rf2(:,6,1)), sqrt(((theta_rf2(:,6,1)-ones(R,1)*delta0)'*(theta_rf2(:,6,1)-ones(R,1)*delta0))/R),        mean(theta_rf2_sd(:,6,1)), mean(theta_rf2_cp(:,6,1)));

fprintf('\n Table 5B \t Bias corrected estimates, number of factors in the y equation is over-specified by 2 \n');
fprintf('\t beta_z      \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',beta_z0,    mean(theta_rf2(:,1,2))-beta_z0,      std(theta_rf2(:,1,2)), sqrt(((theta_rf2(:,1,2)-ones(R,1)*beta_z0)'*(theta_rf2(:,1,2)-ones(R,1)*beta_z0))/R),      mean(theta_rf2_sd(:,1,2)), mean(theta_rf2_cp(:,1,2)));
fprintf('\t beta_y      \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',beta_y0,    mean(theta_rf2(:,2,2))-beta_y0,      std(theta_rf2(:,2,2)), sqrt(((theta_rf2(:,2,2)-ones(R,1)*beta_y0)'*(theta_rf2(:,2,2)-ones(R,1)*beta_y0))/R),      mean(theta_rf2_sd(:,2,2)), mean(theta_rf2_cp(:,2,2)));
fprintf('\t lambda      \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',lambda0,    mean(theta_rf2(:,3,2))-lambda0,      std(theta_rf2(:,3,2)), sqrt(((theta_rf2(:,3,2)-ones(R,1)*lambda0)'*(theta_rf2(:,3,2)-ones(R,1)*lambda0))/R),      mean(theta_rf2_sd(:,3,2)), mean(theta_rf2_cp(:,3,2)));
fprintf('\t alpha_xi    \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',alpha_xi0,  mean(theta_rf2(:,4,2))-alpha_xi0,    std(theta_rf2(:,4,2)), sqrt(((theta_rf2(:,4,2)-ones(R,1)*alpha_xi0)'*(theta_rf2(:,4,2)-ones(R,1)*alpha_xi0))/R),  mean(theta_rf2_sd(:,4,2)), mean(theta_rf2_cp(:,4,2)));
fprintf('\t alpha       \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',alpha0,     mean(theta_rf2(:,5,2))-alpha0,       std(theta_rf2(:,5,2)), sqrt(((theta_rf2(:,5,2)-ones(R,1)*alpha0)'*(theta_rf2(:,5,2)-ones(R,1)*alpha0))/R),        mean(theta_rf2_sd(:,5,2)), mean(theta_rf2_cp(:,5,2)));
fprintf('\t delta       \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',delta0,     mean(theta_rf2(:,6,2))-delta0,       std(theta_rf2(:,6,2)), sqrt(((theta_rf2(:,6,2)-ones(R,1)*delta0)'*(theta_rf2(:,6,2)-ones(R,1)*delta0))/R),        mean(theta_rf2_sd(:,6,2)), mean(theta_rf2_cp(:,6,2)));

fprintf('\n Table 6A \t QML estimates, number of factors in the y equation is over-specified by 3 \n');
fprintf('\t beta_z      \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',beta_z0,    mean(theta_rf3(:,1,1))-beta_z0,      std(theta_rf3(:,1,1)), sqrt(((theta_rf3(:,1,1)-ones(R,1)*beta_z0)'*(theta_rf3(:,1,1)-ones(R,1)*beta_z0))/R),      mean(theta_rf3_sd(:,1,1)), mean(theta_rf3_cp(:,1,1)));
fprintf('\t beta_y      \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',beta_y0,    mean(theta_rf3(:,2,1))-beta_y0,      std(theta_rf3(:,2,1)), sqrt(((theta_rf3(:,2,1)-ones(R,1)*beta_y0)'*(theta_rf3(:,2,1)-ones(R,1)*beta_y0))/R),      mean(theta_rf3_sd(:,2,1)), mean(theta_rf3_cp(:,2,1)));
fprintf('\t lambda      \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',lambda0,    mean(theta_rf3(:,3,1))-lambda0,      std(theta_rf3(:,3,1)), sqrt(((theta_rf3(:,3,1)-ones(R,1)*lambda0)'*(theta_rf3(:,3,1)-ones(R,1)*lambda0))/R),      mean(theta_rf3_sd(:,3,1)), mean(theta_rf3_cp(:,3,1)));
fprintf('\t alpha_xi    \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',alpha_xi0,  mean(theta_rf3(:,4,1))-alpha_xi0,    std(theta_rf3(:,4,1)), sqrt(((theta_rf3(:,4,1)-ones(R,1)*alpha_xi0)'*(theta_rf3(:,4,1)-ones(R,1)*alpha_xi0))/R),  mean(theta_rf3_sd(:,4,1)), mean(theta_rf3_cp(:,4,1)));
fprintf('\t alpha       \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',alpha0,     mean(theta_rf3(:,5,1))-alpha0,       std(theta_rf3(:,5,1)), sqrt(((theta_rf3(:,5,1)-ones(R,1)*alpha0)'*(theta_rf3(:,5,1)-ones(R,1)*alpha0))/R),        mean(theta_rf3_sd(:,5,1)), mean(theta_rf3_cp(:,5,1)));
fprintf('\t delta       \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',delta0,     mean(theta_rf3(:,6,1))-delta0,       std(theta_rf3(:,6,1)), sqrt(((theta_rf3(:,6,1)-ones(R,1)*delta0)'*(theta_rf3(:,6,1)-ones(R,1)*delta0))/R),        mean(theta_rf3_sd(:,6,1)), mean(theta_rf3_cp(:,6,1)));

fprintf('\n Table 6B \t Bias corrected estimates, number of factors in the y equation is over-specified by 3 \n');
fprintf('\t beta_z      \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',beta_z0,    mean(theta_rf3(:,1,2))-beta_z0,      std(theta_rf3(:,1,2)), sqrt(((theta_rf3(:,1,2)-ones(R,1)*beta_z0)'*(theta_rf3(:,1,2)-ones(R,1)*beta_z0))/R),      mean(theta_rf3_sd(:,1,2)), mean(theta_rf3_cp(:,1,2)));
fprintf('\t beta_y      \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',beta_y0,    mean(theta_rf3(:,2,2))-beta_y0,      std(theta_rf3(:,2,2)), sqrt(((theta_rf3(:,2,2)-ones(R,1)*beta_y0)'*(theta_rf3(:,2,2)-ones(R,1)*beta_y0))/R),      mean(theta_rf3_sd(:,2,2)), mean(theta_rf3_cp(:,2,2)));
fprintf('\t lambda      \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',lambda0,    mean(theta_rf3(:,3,2))-lambda0,      std(theta_rf3(:,3,2)), sqrt(((theta_rf3(:,3,2)-ones(R,1)*lambda0)'*(theta_rf3(:,3,2)-ones(R,1)*lambda0))/R),      mean(theta_rf3_sd(:,3,2)), mean(theta_rf3_cp(:,3,2)));
fprintf('\t alpha_xi    \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',alpha_xi0,  mean(theta_rf3(:,4,2))-alpha_xi0,    std(theta_rf3(:,4,2)), sqrt(((theta_rf3(:,4,2)-ones(R,1)*alpha_xi0)'*(theta_rf3(:,4,2)-ones(R,1)*alpha_xi0))/R),  mean(theta_rf3_sd(:,4,2)), mean(theta_rf3_cp(:,4,2)));
fprintf('\t alpha       \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',alpha0,     mean(theta_rf3(:,5,2))-alpha0,       std(theta_rf3(:,5,2)), sqrt(((theta_rf3(:,5,2)-ones(R,1)*alpha0)'*(theta_rf3(:,5,2)-ones(R,1)*alpha0))/R),        mean(theta_rf3_sd(:,5,2)), mean(theta_rf3_cp(:,5,2)));
fprintf('\t delta       \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \t T-SD = %.5f \t CP = %.5f \n',delta0,     mean(theta_rf3(:,6,2))-delta0,       std(theta_rf3(:,6,2)), sqrt(((theta_rf3(:,6,2)-ones(R,1)*delta0)'*(theta_rf3(:,6,2)-ones(R,1)*delta0))/R),        mean(theta_rf3_sd(:,6,2)), mean(theta_rf3_cp(:,6,2)));

fprintf('\n Table 7 \t Estimation assuming exogenous spatial weights matrices \n');
fprintf('\t beta_y      \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \n',beta_y0,    mean(theta_ex(:,1))-beta_y0,      std(theta_ex(:,1)), sqrt(((theta_ex(:,1)-ones(R,1)*beta_y0)'*(theta_ex(:,1)-ones(R,1)*beta_y0))/R));
fprintf('\t lambda      \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \n',lambda0,    mean(theta_ex(:,2))-lambda0,      std(theta_ex(:,2)), sqrt(((theta_ex(:,2)-ones(R,1)*lambda0)'*(theta_ex(:,2)-ones(R,1)*lambda0))/R));
fprintf('\t alpha_xi    \t\t true = %.5f \t bias = %.5f \t E-SD = %.5f \t RMSE = %.5f \n',alpha_xi0,  mean(theta_ex(:,3))-alpha_xi0,    std(theta_ex(:,3)), sqrt(((theta_ex(:,3)-ones(R,1)*alpha_xi0)'*(theta_ex(:,3)-ones(R,1)*alpha_xi0))/R));

save(data_name,'W','X','Y','Z','Gamma_z','Gamma_y','F_z','F_y','factor_error_ER','factor_error_GR','theta','theta_sd','theta_cp','theta_rf1','theta_rf1_sd','theta_rf1_cp','theta_rf2','theta_rf2_sd','theta_rf2_cp','theta_rf3','theta_rf3_sd','theta_rf3_cp','theta_ex')

toc
diary off