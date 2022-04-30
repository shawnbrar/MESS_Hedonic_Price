% This function computes some statistics given estimates
%   input:  theta: parameter vector
%           opt = 1: original parameter
%           opt = 0: the transformed parameter
%
%   output: bc: bias corrected estimator
%           sd: estimated theoretical standard deviation given theta
%               obtained from diagonal elements of C^{-1}/nT.
%
%   10/15/2015

function [bc,sd,gamma_y,f_y,C,SIG] = bias_correction_exog(theta,n,T,Ry,W,X,Y,info,biasopt)

fieldsinf = fieldnames(info);
    nf = length(fieldsinf);
    for i=1:nf
        if strcmp(fieldsinf{i},'boundtyp')
            boundtyp = info.boundtyp; 
        elseif strcmp(fieldsinf{i},'optimtyp')
            optimtyp = info.optimtyp;
        end
    end
%global n T Rz Ry W X Y Z

ky = size(X,2)/T;
beta_y      = theta(1:ky,1);

if strcmp(optimtyp,'Shi_con')==1
    lambda      = theta(ky+1);
    alpha_xi    = theta(ky+2);

    if boundtyp == 1
        sigma_xi = sqrt(alpha_xi);
    elseif boundtyp == 2
        sigma_xi  = 1/alpha_xi;
    end
elseif strcmp(optimtyp,'Shi_unc')==1
    if biasopt == 0
        lambda  = (exp(2*theta(ky+1))-1)/(exp(2*theta(ky+1))+1);  % lambda < 1, > -1
        alpha_xi= exp(theta(ky+2));    % sigma_xi^{-1}

        sigma_xi  = 1/alpha_xi;
    elseif biasopt == 1
        lambda  = theta(ky+1);  
        alpha_xi= theta(ky+2);    
     
        sigma_xi  = 1/alpha_xi;
    end
end

est = [beta_y;lambda;alpha_xi];

%% Matrix G_tilde (nT x nT)
GT = zeros(n*T,n*T);
for t = 1:T
    n0 = 1+(t-1)*n;
    n1 = t*n;
    GT(n0:n1,n0:n1) = W(:,:,t)/(eye(n)-lambda*W(:,:,t));
end

%% Spatially generated regressors
X2 = zeros(n,T);
for t = 1:T
    X2(:,t) = W(:,:,t)/(eye(n)-lambda*W(:,:,t))*(Y(:,t)-lambda*W(:,:,t)*Y(:,t));
end

%% Gamma_y
G = zeros(n,T);
cX_Y = X*kron(beta_y,eye(T));

for t = 1:T
    S = eye(n)-lambda*W(:,:,t);
    G(:,t) = (1/sigma_xi)*(S*Y(:,t)-cX_Y(:,t));
end
[U,D,~] = svd(G);
D = diag(D).*diag(D);
[~,I] = sort(D);
Gamma_y = zeros(n,Ry);

for c = 1:Ry
    Gamma_y(:,c) = U(:,I(end-c+1));
end
P_Gy = Gamma_y*Gamma_y';
M_Gy = eye(n) - P_Gy;

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

%% Bias in C^(1) (6x1)
phi = [
        zeros(ky,1);
        -(1/sqrt(n*T))*trace(kron(P_Fy,eye(n))*GT+kron(eye(T),P_Gy)*GT);
        (sqrt(n/T)+sqrt(T/n))*Ry*sigma_xi;
        ];

%% Matrix C
% SIGMA matrix
SIG = zeros(ky+2); % kz+ky+2+J+p

% Elements corresponding to \beta_y\beta_y
for k1 = 1:ky
    X_y_k1 = X(:,T*(k1-1)+1:k1*T);
    for k2 = 1:ky
        X_y_k2 = X(:,T*(k2-1)+1:k2*T);
        SIG(k1,k2) = trace(M_Fy*X_y_k1'*M_Gy*X_y_k2)/(n*T*sigma_xi^2);
    end
end

% Elements corresponding to \beta_y\lambda
for k1 = 1:ky
    X_y_k1 = X(:,T*(k1-1)+1:k1*T);
    SIG(ky+1,k1) = trace(M_Fy*X_y_k1'*M_Gy*X2)/(n*T*sigma_xi^2);
    SIG(k1,ky+1) = trace(M_Fy*X_y_k1'*M_Gy*X2)/(n*T*sigma_xi^2);
end

% Elements corresponding to \lambda\lambda
SIG(ky+1,ky+1) = trace(M_Fy*X2'*M_Gy*X2)/(n*T*sigma_xi^2) + ...
                                                       trace(GT*GT)/(n*T);
                                                   
% Elements corresponding to \sigma_xi\beta_y and \alpha\beta_y
for k1 = 1:ky
    SIG(ky+2,k1) = 0;
    SIG(k1,ky+2) = 0;
end

% Elements corresponding to \sigma_xi\lambada and \alpha\lambda
SIG(ky+1,ky+2) = -trace(GT)*2*sigma_xi/(n*T);
SIG(ky+2,ky+1) = -trace(GT)*2*sigma_xi/(n*T);


% Elements corresponding to \sigma_xi\alpha, \sigma_xi\sigma_xi and \alpha\alpha
SIG(ky+2,ky+2) = 2*sigma_xi^2;

% Case of normally distributed error
C = SIG;
% C = zeros(6,6);
% 
% C(1,1) = trace(M_Fz*X'*M_Gz*X)/(n*T*sigma_e^2) + trace(M_Fy*X'*delta*M_Gy*delta*X)/(n*T*sigma_xi^2);
% C(1,2) = -trace(M_Fy*X'*delta*M_Gy*X)/(n*T*sigma_xi^2);
% C(1,3) = -trace(M_Fy*X'*delta*M_Gy*X2)/(n*T*sigma_xi^2);
% C(1,4) = 0;
% C(1,5) = 0;
% C(1,6) = -trace(M_Gy*delta*X*M_Fy*(Z-X_Z*kron(beta_z,eye(T)))')/(n*T*sigma_xi^2);
% 
% C(2,1) = C(1,2);
% C(2,2) = trace(M_Fy*X'*M_Gy*X)/(n*T*sigma_xi^2);
% C(2,3) = trace(M_Fy*X'*M_Gy*X2)/(n*T*sigma_xi^2);
% C(2,4) = 0;
% C(2,5) = 0;
% C(2,6) = -trace(M_Gy*X*M_Fy*(Z-X_Z*kron(beta_z,eye(T)))')/(n*T*sigma_xi^2);
% 
% C(3,1) = C(1,3);
% C(3,2) = C(2,3);
% C(3,3) = trace(M_Fy*X2'*M_Gy*X2)/(n*T*sigma_xi^2) + trace(GT*GT)/(n*T);
% C(3,4) = -trace(GT)*2*sigma_xi/(n*T);
% C(3,5) = 0;
% C(3,6) = -trace(M_Gy*X2*M_Fy*(Z-X_Z*kron(beta_z,eye(T)))')/(n*T*sigma_xi^2);
% 
% C(4,1) = C(1,4);
% C(4,2) = C(2,4);
% C(4,3) = C(3,4);
% C(4,4) = 2*sigma_xi^2;
% C(4,5) = 0;
% C(4,6) = 0;
% 
% C(5,1) = C(1,5);
% C(5,2) = C(2,5);
% C(5,3) = C(3,5);
% C(5,4) = C(4,5);
% C(5,5) = 2*sigma_e^2;
% C(5,6) = 0;
% 
% C(6,1) = C(1,6);
% C(6,2) = C(2,6);
% C(6,3) = C(3,6);
% C(6,4) = C(4,6);
% C(6,5) = C(5,6);
% C(6,6) = trace(M_Gy*(Z-X_Z*kron(beta_z,eye(T)))*M_Fy*(Z-X_Z*kron(beta_z,eye(T)))')/(n*T*sigma_xi^2);
% 
%% bias corrected estimator
%size(phi)
%size(C)

bc = est - (C\phi)/sqrt(n*T);

% Theoretical standard deviation for the original estimator
%  obtained from diagonal elements of C^{-1}/nT.
CT = invpd(C)/(n*T);
sd = sqrt(diag(CT));

% sd = zeros(6,1);
% sd(1) = sqrt(CT(1,1));
% sd(2) = sqrt(CT(2,2));
% sd(3) = sqrt(CT(3,3));
% sd(4) = sqrt(CT(4,4));
% sd(5) = sqrt(CT(5,5));
% sd(6) = sqrt(CT(6,6));

 % factor loadings
    gamma_y = Gamma_y;
    
    % time factors, obtained from cross sectional regressions
    f_y = zeros(T,Ry);
    for t = 1:T
        res = (eye(n)-lambda*squeeze(W(:,:,t)))*Y(:,t)- cX_Y(:,t);
        f_y(t,:) = ((gamma_y'*gamma_y)\(gamma_y'*res))'; 
    end

end
