%==========================================================================
% Spatial dynamic panel data models with interactive fixed effects and 
%   endogenous weights matrices
% Last updated: 2017/6/28
% This program implements the models in Shi and Lee (2017) "A spatial panel
%   data model with time varying endogenous weights matrices and common
%   factors", Regional Science and Urban Economics, 2017,
%   10.1016/j.regsciurbeco.2017.03.007
%
% This code is offered with no guarantee.
% I appreciate for any feedback on issues.
% 
%==========================================================================
% Author
%   Wei Shi, Institute for Economic and Social Research, Jinan University
%   Email: shiweiemail@gmail.com
%
%==========================================================================
% Method
% The profile objective function (Eq.6) is maximized using fminunc.
%
%==========================================================================
% Input parameters
%   Y: NxT matrix of outcomes
%   W: NxNxT time varying spatial weights matrix in the outcome equation
%       W(:,:,t) = W(t)
%   Xy: KyxNxT regressors in y
%   Xz: KzxNxT regressors in z
%   Z: NxT regressors that construct the spatial weights (see section
%   2.2)
%   Ry: number of factors in y
%   Rz: number of factors in z 
%   TolX: stoping criteria on beta
%
% Output parameters
%   beta: QML estimate with entries corresponding to (see Eq.1,3 for model
%       specification)
%       beta(1:Kz):         Xz
%       beta(Kz+1:Ky+Kz):   Xy 
%       beta(Ky+Kz+1):      WY
%       beta(Ky+Kz+2):      sigma_xi^{-1}
%       beta(Ky+Kz+3):      sigma_epsilon^{-1} (see Assumption 7)
%       beta(Ky+Kz+4):      delta
%   bcorr: bias correction term (see Section 4.3 of the paper)
%       the bias corrected estimator is beta-bcorr
%   Vbeta: estimated variance-covariance matrix (assuming normal errors)
%   gamma_z: factor loadings in the z equation
%   f_z: time factors in the z equation
%   gamma_y: factor loadings in the y equation
%   f_y: time factors in the y equation
%   exitflag: Positive exit flags correspond to successful outcomes.

function [beta,bcorr,Vbeta,gamma_z,f_z,gamma_y,f_y,exitflag] = ESFactors(Y,W,Xy,Xz,Z,Ry,Rz,TolX)
    N = size(Y,1);
    T = size(Y,2);
    Ky = size(Xy,1);
    Kz = size(Xz,1);
    
    opt = optimoptions('fminunc','TolX',TolX,'MaxFunEvals',100000,'MaxIter',100000,'Display','notify','Algorithm','quasi-newton');
    % try multiple starting values
    fval0 = 1/eps;
    for i = 1:9
        start = randn(Kz+Ky+4,1);
        [est,fval,exitflag0] = fminunc(@Q_obj,start,opt,Y,W,Xy,Xz,Z,Ry,Rz);
        if fval < fval0 && exitflag0 > 0
            fval0 = fval;
            beta = est;
            exitflag = exitflag0;
        end
    end
    
    beta_z  = beta(1:Kz);
    beta_y  = beta(Kz+1:Ky+Kz);
    lambda  = (exp(2*beta(Ky+Kz+1))-1)/(exp(2*beta(Ky+Kz+1))+1);  % lambda < 1, > -1
    alpha_xi= exp(beta(Ky+Kz+2));    % sigma_xi^{-1}
    alpha   = exp(beta(Ky+Kz+3));    % sigma_epsilon^{-1}
    delta   = beta(Ky+Kz+4);
    sigma_e     = 1/alpha;
    sigma_xi    = 1/alpha_xi;    
    
    Xzsum = zeros(N,T);
    for k = 1:Kz
        Xzsum = Xzsum + beta_z(k)*squeeze(Xz(k,:,:));
    end
    Xysum = zeros(N,T);
    for k = 1:Ky
        Xysum = Xysum + beta_y(k)*squeeze(Xy(k,:,:));
    end          
    
    %% Matrix G_tilde (NT x NT)
    GT = zeros(N*T,N*T);
    for t = 1:T
        n0 = 1+(t-1)*N;
        n1 = t*N;
        GT(n0:n1,n0:n1) = W(:,:,t)/(eye(N)-lambda*W(:,:,t));
    end
    
    %% Spatially generated regressors
    X_SAR = zeros(N,T);
    for t = 1:T
        X_SAR(:,t) = W(:,:,t)*Y(:,t);
    end
    
    %% Projection matrices
    %% Gamma_z
    M = (1/sigma_e)*(Z-Xzsum);
    [U,D,~] = svd(M);
    D = diag(D).*diag(D);
    [~,I] = sort(D);
    Gamma_z = zeros(N,Rz);
    for c = 1:Rz
        Gamma_z(:,c) = U(:,I(end-c+1));
    end
    P_Gz = Gamma_z*Gamma_z';
    M_Gz = eye(N) - P_Gz;

    %% Time factors F_z
    [U,D,~] = svd(M');
    D = diag(D).*diag(D);
    [~,I] = sort(D);
    F_z = zeros(T,Rz);
    for c = 1:Rz
        F_z(:,c) = U(:,I(end-c+1));
    end
    P_Fz = F_z*F_z';
    M_Fz = eye(T) - P_Fz;

    %% Gamma_y
    G = zeros(N,T);
    for t = 1:T
        S = eye(N)-lambda*squeeze(W(:,:,t));
        G(:,t) = (1/sigma_xi)*(S*Y(:,t)-Xysum(:,t)-(Z(:,t)-Xzsum(:,t))*delta);
    end
    [U,D,~] = svd(G);
    D = diag(D).*diag(D);
    [~,I] = sort(D);
    Gamma_y = zeros(N,Ry);
    for c = 1:Ry
        Gamma_y(:,c) = U(:,I(end-c+1));
    end
    P_Gy = Gamma_y*Gamma_y';
    M_Gy = eye(N) - P_Gy;

    %% Time factors F_y
    [U,D,~] = svd(G');
    D = diag(D).*diag(D);
    [~,I] = sort(D);
    F_y = zeros(T,Ry);
    for c = 1:Ry
        F_y(:,c) = U(:,I(end-c+1));
    end
    P_Fy = F_y*F_y';
    M_Fy = eye(T) - P_Fy;    
    
    %% Matrix C
    C = zeros(Ky+Kz+4,Ky+Kz+4);

    for i = 1:Kz
        for j = 1:Kz
            C(i,j) = trace(M_Fz*squeeze(Xz(i,:,:))'*M_Gz*squeeze(Xz(j,:,:)))/(N*T*sigma_e^2) + trace(M_Fy*squeeze(Xz(i,:,:))'*delta*M_Gy*delta*squeeze(Xz(j,:,:)))/(N*T*sigma_xi^2);
        end
        for j = 1:Ky
            C(i,Kz+j) = -trace(M_Fy*squeeze(Xz(i,:,:))'*delta*M_Gy*squeeze(Xy(j,:,:)))/(N*T*sigma_xi^2);
        end
        C(i,Ky+Kz+1) = -trace(M_Fy*squeeze(Xz(i,:,:))'*delta*M_Gy*X_SAR)/(N*T*sigma_xi^2);
        C(i,Ky+Kz+2) = 0;
        C(i,Ky+Kz+3) = 0;
        C(i,Ky+Kz+4) = -trace(M_Gy*delta*squeeze(Xz(i,:,:))*M_Fy*(Z-Xzsum)')/(N*T*sigma_xi^2);
    end
    
    for i = 1:Ky
        for j = 1:Kz
            C(Kz+i,j) = C(j,Kz+i);
        end
        for j = 1:Ky
            C(Kz+i,Kz+j) = trace(M_Fy*squeeze(Xy(i,:,:))'*M_Gy*squeeze(Xy(j,:,:)))/(N*T*sigma_xi^2);
        end
        C(Kz+i,Ky+Kz+1) = trace(M_Fy*squeeze(Xy(i,:,:))'*M_Gy*X_SAR)/(N*T*sigma_xi^2);
        C(Kz+i,Ky+Kz+2) = 0;
        C(Kz+i,Ky+Kz+3) = 0;
        C(Kz+i,Ky+Kz+4) = -trace(M_Gy*squeeze(Xy(i,:,:))*M_Fy*(Z-Xzsum)')/(N*T*sigma_xi^2);
    end
    
    % lambda
    for j = 1:Kz
        C(Kz+Ky+1,j) = C(j,Kz+Ky+1);
    end
    for j = 1:Ky
        C(Kz+Ky+1,Kz+j) = C(Kz+j,Kz+Ky+1);
    end
    C(Kz+Ky+1,Ky+Kz+1) = trace(M_Fy*X_SAR'*M_Gy*X_SAR)/(N*T*sigma_xi^2) + trace(GT*GT)/(N*T);
    C(Kz+Ky+1,Ky+Kz+2) = -trace(GT)*2*sigma_xi/(N*T);
    C(Kz+Ky+1,Ky+Kz+3) = 0;
    C(Kz+Ky+1,Ky+Kz+4) = -trace(M_Gy*X_SAR*M_Fy*(Z-Xzsum)')/(N*T*sigma_xi^2);

    for j = 1:Kz
        C(Kz+Ky+2,j) = C(j,Kz+Ky+2);
    end
    for j = 1:Ky
        C(Kz+Ky+2,Kz+j) = C(Kz+j,Kz+Ky+2);
    end
    C(Kz+Ky+2,Ky+Kz+1) = C(Ky+Kz+1,Kz+Ky+2);
    C(Kz+Ky+2,Ky+Kz+2) = 2*sigma_xi^2;
    C(Kz+Ky+2,Ky+Kz+3) = 0;
    C(Kz+Ky+2,Ky+Kz+4) = 0;
    
    for j = 1:Kz
        C(Kz+Ky+3,j) = C(j,Kz+Ky+3);
    end
    for j = 1:Ky
        C(Kz+Ky+3,Kz+j) = C(Kz+j,Kz+Ky+3);
    end
    C(Kz+Ky+3,Ky+Kz+1) = C(Ky+Kz+1,Kz+Ky+3);
    C(Kz+Ky+3,Ky+Kz+2) = C(Ky+Kz+2,Kz+Ky+3);
    C(Kz+Ky+3,Ky+Kz+3) = 2*sigma_e^2;
    C(Kz+Ky+3,Ky+Kz+4) = 0;

    for j = 1:Kz
        C(Kz+Ky+4,j) = C(j,Kz+Ky+4);
    end
    for j = 1:Ky
        C(Kz+Ky+4,Kz+j) = C(Kz+j,Kz+Ky+4);
    end
    C(Kz+Ky+4,Ky+Kz+1) = C(Ky+Kz+1,Kz+Ky+4);
    C(Kz+Ky+4,Ky+Kz+2) = C(Ky+Kz+2,Kz+Ky+4);
    C(Kz+Ky+4,Ky+Kz+3) = C(Kz+Ky+3,Ky+Kz+4);
    C(Kz+Ky+4,Ky+Kz+4) = trace(M_Gy*(Z-Xzsum)*M_Fy*(Z-Xzsum)')/(N*T*sigma_xi^2);

    %% Bias in C^(1) ((Ky+Kz+4)x1)
    phi = [
            zeros(Kz+Ky,1);
            -(1/sqrt(N*T))*trace(kron(P_Fy,eye(N))*GT+kron(eye(T),P_Gy)*GT);
            (sqrt(N/T)+sqrt(T/N))*Ry*sigma_xi;
            (sqrt(N/T)+sqrt(T/N))*Rz*sigma_e;
            0
            ];    
    
    % parameter vector: gamma;rho;beta;lambda;alpha
    beta = [beta_z;beta_y;lambda;alpha_xi;alpha;delta];
    bcorr = (C\phi)/sqrt(N*T);
    Vbeta = inv(C)/(N*T);
    
    % factor loadings
    gamma_z = Gamma_z;
    gamma_y = Gamma_y;
    
    % time factors, obtained from cross sectional regressions
    f_y = zeros(T,Ry);
    for t = 1:T
        res = (eye(N)-lambda*squeeze(W(:,:,t)))*Y(:,t)-Xysum(:,t)-(Z(:,t)-Xzsum(:,t))*delta;
        f_y(t,:) = ((gamma_y'*gamma_y)\(gamma_y'*res))'; 
    end
    f_z = ((gamma_z'*gamma_z)\(gamma_z'*(Z-Xzsum)))'; 
return

function f = Q_obj(theta,Y,W,Xy,Xz,Z,Ry,Rz)
    % This function evaluates the likelihood function Q_nT(theta), see
    % Eq.6.
    % Ry(Rz) factors are concentrated out from y(z).
    N = size(Y,1);
    T = size(Y,2);
    Ky = size(Xy,1);
    Kz = size(Xz,1);
   
    beta_z  = theta(1:Kz);
    beta_y  = theta(Kz+1:Ky+Kz);
    lambda  = (exp(2*theta(Ky+Kz+1))-1)/(exp(2*theta(Ky+Kz+1))+1);  % lambda < 1, > -1
    alpha_xi= exp(theta(Ky+Kz+2));    % sigma_xi^{-1}
    alpha   = exp(theta(Ky+Kz+3));    % sigma_epsilon^{-1}
    delta   = theta(Ky+Kz+4);

    sigma_eps_sq = 1/alpha^2;
    sigma_xi_sq  = 1/alpha_xi^2;

    Xzsum = zeros(N,T);
    for k = 1:Kz
        Xzsum = Xzsum + beta_z(k)*squeeze(Xz(k,:,:));
    end
    Xysum = zeros(N,T);
    for k = 1:Ky
        Xysum = Xysum + beta_y(k)*squeeze(Xy(k,:,:));
    end      
    
    ll = -(1/2)*log(sigma_eps_sq)-(1/2)*log(sigma_xi_sq);
    for t = 1:T 
        S = eye(N)-lambda*squeeze(W(:,:,t));
        ll = ll + log(det(S))/(N*T);
    end      
    
    D = (1/sqrt(sigma_eps_sq))*(Z-Xzsum);
    G = zeros(N,T);
    for t = 1:T
        S = eye(N)-lambda*squeeze(W(:,:,t));
        G(:,t) = (1/sqrt(sigma_xi_sq))*(S*Y(:,t)-Xysum(:,t)-(Z(:,t)-Xzsum(:,t))*delta);
    end    
    
if sum(sum(isnan(D))) ~= 0 || sum(sum(isinf(D))) ~= 0 || sum(sum(isnan(G))) ~= 0 || sum(sum(isinf(G))) ~= 0 || isnan(ll) == 1 || isinf(ll) == 1
    f = 4.503e+15;
else
    [~,D1,~] = svd(D);
    D1 = diag(D1).*diag(D1);
    D1 = sort(D1);

    [~,G1,~] = svd(G);
    G1 = diag(G1).*diag(G1);
    G1 = sort(G1);

    ll = ll-(1/2)*sum(D1(1:end-Rz))/(N*T)-(1/2)*sum(G1(1:end-Ry))/(N*T);
    f = -ll;
end    

return