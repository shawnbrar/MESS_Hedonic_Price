% Monte Carlo experiment for "ESFactors.m"

tic

%% Set Monte Carlo design parameters
n   = 25;
T   = 25;
Rz0 = 1;      % true number of factors
Ry0 = 2;      % true number of factors in the y equation
SNR = 1;      % variance of v and e, 1,2,4,9
R   = 10;   % number of Monte Carlo iterations
rng(20151014);

beta_z0  = [-1;1];
beta_y0  = [1;-1];
lambda0  = 0.2;
sigma_v0 = sqrt(16/12*SNR);
sigma_e0 = sqrt(16/12*SNR);
rho0     = 0.2;

%% Parameters
sigma_xi0 = sqrt((1-rho0^2)*sigma_v0^2);
alpha_xi0 = 1/sigma_xi0;
alpha0 = 1/sigma_e0;
delta0 = rho0*sigma_v0/sigma_e0;
m      = min(n,T);

%% Generate Data
NT      = n*T;
Gamma_y = rand(n,Ry0)*4-ones(n,Ry0)*2;
Gamma_z = rand(n,Rz0)*4-ones(n,Rz0)*2;
F_y     = rand(T,Ry0)*4-ones(T,Ry0)*2;
F_z     = F_y(:,1:Rz0);

Xz = zeros(2,n,T);
Xz(1,:,:) = rand(n,T)*4-ones(n,T)*2 + (Gamma_z*F_z' + Gamma_z*ones(Rz0,1)*ones(1,T) + ones(n,1)*(F_z*ones(Rz0,1))')/3;
Xz(2,:,:) = rand(n,T)*4-ones(n,T)*2;
Xy = zeros(2,n,T);
Xy(1,:,:) = rand(n,T)*4-ones(n,T)*2 + (Gamma_y*F_y' + Gamma_y*ones(Ry0,1)*ones(1,T) + ones(n,1)*(F_y*ones(Ry0,1))')/3;
Xy(2,:,:) = rand(n,T)*4-ones(n,T)*2;
mu      = [0,0];
Sigma   = [sigma_v0^2,rho0*sigma_v0*sigma_e0;rho0*sigma_v0*sigma_e0,sigma_e0^2];

%% Create the spatial weights matrix from rook adjacency matrix
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

%% Monte Carlo Iterations
theta_e     = zeros(8,R);
theta_bc    = zeros(8,R);
cp          = zeros(8,R);

for mc_r = 1:R
    V   = zeros(n,T);
    E   = zeros(n,T);
    Y   = zeros(n,T);
    W   = zeros(n,n,T);
    
    Error = mvnrnd(mu,Sigma,NT);
    for t = 1:T
        V(:,t) = Error(1+(t-1)*n:t*n,1);
        E(:,t) = Error(1+(t-1)*n:t*n,2);
    end
    Z = squeeze(Xz(1,:,:))*beta_z0(1) + squeeze(Xz(2,:,:))*beta_z0(2) + Gamma_z*F_z' + E;    
    
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
        Y(:,t) = S0\(Xy(1,:,t)'*beta_y0(1) + Xy(2,:,t)'*beta_y0(2) + Gamma_y*F_y(t,:)' + V(:,t));
    end
    
    opt = optimoptions('fminunc','TolFun',1e-9,'TolX',1e-9,'MaxFunEvals',100000,'MaxIter',100000,'Display','notify','Algorithm','quasi-newton');
    
    fprintf('Monte Carlo Iterations \t round %.0f \t %.2f%% \t \n',mc_r,(mc_r/R)*100);
    [beta,bcorr,Vbeta,gamma_z,f_z,gamma_y,f_y,exitflag] = ESFactors(Y,W,Xy,Xz,Z,Ry0,Rz0,1e-9);
    theta_e(:,mc_r) = beta;
    theta_bc(:,mc_r)= beta-bcorr;    
    
    % coverage probability
    if (1-normcdf(abs((beta_z0(1)-theta_bc(1,mc_r))/sqrt(Vbeta(1,1)))))*2 >= 0.05
        cp(1,mc_r) = 1;
    end
    if (1-normcdf(abs((beta_z0(2)-theta_bc(2,mc_r))/sqrt(Vbeta(2,2)))))*2 >= 0.05
        cp(2,mc_r) = 1;
    end
    if (1-normcdf(abs((beta_y0(1)-theta_bc(3,mc_r))/sqrt(Vbeta(3,3)))))*2 >= 0.05
        cp(3,mc_r) = 1;
    end
    if (1-normcdf(abs((beta_y0(2)-theta_bc(4,mc_r))/sqrt(Vbeta(4,4)))))*2 >= 0.05
        cp(4,mc_r) = 1;
    end    
    if (1-normcdf(abs((lambda0-theta_bc(5,mc_r))/sqrt(Vbeta(5,5)))))*2 >= 0.05
        cp(5,mc_r) = 1;
    end
    if (1-normcdf(abs((alpha_xi0-theta_bc(6,mc_r))/sqrt(Vbeta(6,6)))))*2 >= 0.05
        cp(6,mc_r) = 1;
    end
    if (1-normcdf(abs((alpha0-theta_bc(7,mc_r))/sqrt(Vbeta(7,7)))))*2 >= 0.05
        cp(7,mc_r) = 1;
    end
    if (1-normcdf(abs((delta0-theta_bc(8,mc_r))/sqrt(Vbeta(8,8)))))*2 >= 0.05
        cp(8,mc_r) = 1;
    end          
end

fprintf('\n Simulation results: \t n = %.0f T = %.0f \n',n,T);
fprintf('\t beta_z1 \t true = %.5f \t QMLE bias = %.5f \t bias-corrected estimator bias = %.5f \t CP = %.5f\n',beta_z0(1),mean(theta_e(1,:))-beta_z0(1),mean(theta_bc(1,:))-beta_z0(1),mean(cp(1,:)));
fprintf('\t beta_z2 \t true = %.5f \t QMLE bias = %.5f \t bias-corrected estimator bias = %.5f \t CP = %.5f\n',beta_z0(2),mean(theta_e(2,:))-beta_z0(2),mean(theta_bc(2,:))-beta_z0(2),mean(cp(2,:)));
fprintf('\t beta_y1 \t true = %.5f \t QMLE bias = %.5f \t bias-corrected estimator bias = %.5f \t CP = %.5f\n',beta_y0(1),mean(theta_e(3,:))-beta_y0(1),mean(theta_bc(3,:))-beta_y0(1),mean(cp(3,:)));
fprintf('\t beta_y2 \t true = %.5f \t QMLE bias = %.5f \t bias-corrected estimator bias = %.5f \t CP = %.5f\n',beta_y0(2),mean(theta_e(4,:))-beta_y0(2),mean(theta_bc(4,:))-beta_y0(2),mean(cp(4,:)));
fprintf('\t lambda \t true = %.5f \t QMLE bias = %.5f \t bias-corrected estimator bias = %.5f \t CP = %.5f\n',lambda0,mean(theta_e(5,:))-lambda0,mean(theta_bc(5,:))-lambda0,mean(cp(5,:)));
fprintf('\t alpha_xi \t true = %.5f \t QMLE bias = %.5f \t bias-corrected estimator bias = %.5f \t CP = %.5f\n',alpha_xi0,mean(theta_e(6,:))-alpha_xi0,mean(theta_bc(6,:))-alpha_xi0,mean(cp(6,:)));
fprintf('\t alpha \t true = %.5f \t QMLE bias = %.5f \t bias-corrected estimator bias = %.5f \t CP = %.5f\n',alpha0,mean(theta_e(7,:))-alpha0,mean(theta_bc(7,:))-alpha0,mean(cp(7,:)));
fprintf('\t delta \t true = %.5f \t QMLE bias = %.5f \t bias-corrected estimator bias = %.5f \t CP = %.5f\n',delta0,mean(theta_e(8,:))-delta0,mean(theta_bc(8,:))-delta0,mean(cp(8,:)));

toc

function W = row_normalize(W0)

n = length(W0(:,1));
W = zeros(n,n);

W_sum = zeros(n,1);
for i = 1:n
    for j = 1:n
        W_sum(i) = W_sum(i) + W0(i,j);
    end
end

for i = 1:n
    if W_sum(i) ~= 0
        for j = 1:n
            W(i,j) = W0(i,j)/W_sum(i);
        end
    end
end
end