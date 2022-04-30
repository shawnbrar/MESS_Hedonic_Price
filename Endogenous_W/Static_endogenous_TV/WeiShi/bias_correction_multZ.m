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

function [bc,sd,C,SIG] = bias_correction_multZ(b_est,n,T,Rz,Ry,W,X_Y,Y,Z,X_Z,opt)

%b_est = est;
%global n T Rz Ry W X Y Z
kz = (size(X_Z,2)/T)*size(X_Z,3);  %the number of explanatory variables in Z equation
ky = size(X_Y,2)/T;                %the number of explanatory variables in y equation
p  = size(Z,3);                                %the number of dependent variables in the Z equation
J  = p + (p*(p-1))/2;                          %the number of distinct elements in Sigma_epsilon

beta_z     = reshape(b_est(1:kz,1),size(X_Z,2)/T,size(X_Z,3));
beta_y     = b_est(kz+1:kz+ky,1);

if opt == 0
    lambda      = 0.9*(exp(2*b_est(kz+ky+1))-1)/(exp(2*b_est(kz+ky+1))+1);  % lambda < 0.9, > -0.9
    alpha_xi    = exp(b_est(kz+ky+2));    % sigma_xi^{-1}
    alpha       = exp(b_est(kz+ky+3:kz+ky+2+J));    % sigma_epsilon^{-1}
elseif opt == 1
    lambda      = b_est(kz+ky+1);
    alpha_xi    = b_est(kz+ky+2);
    alpha       = b_est(kz+ky+3:kz+ky+2+J);
end

delta       = b_est(kz+ky+2+J+1:end);


%------- See Assumption 7
%sigma_e     = 1/alpha;
sigma_e_inv = zeros(p,p);
lng = p-1:-1:0;
binf = 1;
for p_i = 1:p
    sigma_e_inv(p_i:p,p_i:p) = toeplitz(alpha(binf:binf+lng(p_i),1));%
    binf = binf+lng(p_i) +1;
end
sigma_e_sq = inv(sigma_e_inv*sigma_e_inv);
sigma_e = sqrtm(sigma_e_sq);
%-----------------------------------


%------- See page 7
sigma_xi  = 1/alpha_xi^2;
%--------------------------


b_est = [reshape(beta_z,kz,1);beta_y;lambda;alpha_xi;alpha;delta];

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

%% Projection matrices
%% Gamma_z

for l = 1:p
    cZ(:,:,l) = Z(:,:,l)-X_Z(:,:,l)*kron(beta_z(:,l),eye(T));
end

M = (kron(sigma_e_inv,eye(n)))*[cZ(:,:,1);cZ(:,:,2)];  % no sqrt because here it is not sigma_e_sq ;)
[U,D,~] = svd(M);
D = diag(D).*diag(D);
[~,I] = sort(D);
Gamma_z = zeros(n*p,Rz);
for c = 1:Rz
    Gamma_z(:,c) = U(:,I(end-c+1));   %eigenvectors corresponding to the Rz largest eigenvalues
end
P_Gz = Gamma_z*Gamma_z'; % Normally, P_Gz = Gamma_z*(Gamma_z'*Gamma_z)^-1*Gamma_z'. However as U is orthonormal (Gamma_z'*Gamma_z)=I
M_Gz = eye(n*p) - P_Gz;

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
G = zeros(n,T);
%cZ = Z-X_Z*kron(beta_z,eye(T));
cX_Y = X_Y*kron(beta_y,eye(T));

for t = 1:T
    S = eye(n)-lambda*W(:,:,t);
    G(:,t) = (1/sigma_xi)*(S*Y(:,t)-cX_Y(:,t)-kron(delta',eye(n))*[cZ(:,t,1);cZ(:,t,2)]);
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
        zeros(kz+ky,1);
        -(1/sqrt(n*T))*trace(kron(P_Fy,eye(n))*GT+kron(eye(T),P_Gy)*GT);
        (sqrt(n/T)+sqrt(T/n))*Ry*sigma_xi;
        Rz*sqrt(n/T)*trace(sigma_e*(sigma_e_inv==alpha(1))) + sqrt(n/T)*trace((kron(sigma_e*(sigma_e_inv==alpha(1)),eye(n)))*P_Gz)
        Rz*sqrt(n/T)*trace(sigma_e*(sigma_e_inv==alpha(2))) + sqrt(n/T)*trace((kron(sigma_e*(sigma_e_inv==alpha(2)),eye(n)))*P_Gz)
        Rz*sqrt(n/T)*trace(sigma_e*(sigma_e_inv==alpha(3))) + sqrt(n/T)*trace((kron(sigma_e*(sigma_e_inv==alpha(3)),eye(n)))*P_Gz)
        zeros(p,1)
        ];        %(sqrt(n/T)+sqrt(T/n))*Rz*sigma_e;

%% Matrix C
% SIGMA matrix
SIG = zeros(kz+ky+2+J+p); % kz+ky+2+J+p


% Rewrite X_Z
X_ZT = blockdiag(X_Z(:,:,1),X_Z(:,:,2));
% Elements corresponding to \beta_z\beta_z
for k1 = 1:kz
    X_z_k1 = X_ZT(:,T*(k1-1)+1:k1*T);
    for k2 = 1:kz
        X_z_k2 = X_ZT(:,T*(k2-1)+1:k2*T);
        SIG(k1,k2) = trace(M_Fz*X_z_k1'*kron(inv(sigma_e),eye(n))*M_Gz*kron(inv(sigma_e),eye(n))*X_z_k2)/(n*T) + ...
             trace(M_Fy*X_z_k1'*kron(delta,eye(n))*M_Gy*kron(delta',eye(n))*X_z_k2)/(n*T*sigma_xi^2);
    end
end

% Elements corresponding to \beta_z\beta_y
for k1 = 1:kz
    X_z_k1 = X_ZT(:,T*(k1-1)+1:k1*T);
    for k2 = 1:ky
        X_y_k2 = X_Y(:,T*(k2-1)+1:k2*T);
        SIG(kz+k2,k1) = -trace(M_Fy*X_z_k1'*kron(delta,eye(n))*M_Gy*X_y_k2)/...
                                                        (n*T*sigma_xi^2);
        SIG(k1,kz+k2) = SIG(kz+k2,k1);
    end
end

% Elements corresponding to \beta_y\beta_y
for k1 = 1:ky
    X_y_k1 = X_Y(:,T*(k1-1)+1:k1*T);
    for k2 = 1:ky
        X_y_k2 = X_Y(:,T*(k2-1)+1:k2*T);
        SIG(kz+k1,kz+k2) = trace(M_Fy*X_y_k1'*M_Gy*X_y_k2)/(n*T*sigma_xi^2);
    end
end

% Elements corresponding to \beta_z\lambda
for k1 = 1:kz
    X_z_k1 = X_ZT(:,T*(k1-1)+1:k1*T);
    SIG(kz+ky+1,k1) = -trace(M_Fy*X_z_k1'*(kron(delta',eye(n))')*M_Gy*X2)/(n*T*sigma_xi^2);
    SIG(k1,kz+ky+1) = -trace(M_Fy*X_z_k1'*(kron(delta',eye(n))')*M_Gy*X2)/(n*T*sigma_xi^2);
end

% Elements corresponding to \beta_y\lambda
for k1 = 1:ky
    X_y_k1 = X_Y(:,T*(k1-1)+1:k1*T);
    SIG(kz+ky+1,kz+k1) = trace(M_Fy*X_y_k1'*M_Gy*X2)/(n*T*sigma_xi^2);
    SIG(kz+k1,kz+ky+1) = trace(M_Fy*X_y_k1'*M_Gy*X2)/(n*T*sigma_xi^2);
end

% Elements corresponding to \lambda\lambda
SIG(kz+ky+1,kz+ky+1) = trace(M_Fy*X2'*M_Gy*X2)/(n*T*sigma_xi^2) + ...
                                                       trace(GT*GT)/(n*T);

% Elements corresponding to \sigma_xi\beta_z and \alpha\beta_z
for k1 = 1:kz
    SIG(kz+ky+2,k1) = 0;
    SIG(k1,kz+ky+2) = 0;
    SIG(kz+ky+3:kz+ky+2+J,k1) = 0;  % These two last lines will change in terms of dimension
    SIG(k1,kz+ky+3:kz+ky+2+J) = 0;
end

% Elements corresponding to \sigma_xi\beta_y and \alpha\beta_z
for k1 = 1:ky
    SIG(kz+ky+2,kz+k1) = 0;
    SIG(kz+k1,kz+ky+2) = 0;
    SIG(kz+ky+3:kz+ky+2+J,kz+k1) = 0;  % These two last lines will change in terms of dimension
    SIG(kz+k1,kz+ky+3:kz+ky+2+J) = 0;
end

% Elements corresponding to \sigma_xi\lambda and \alpha\lambda
SIG(kz+ky+1,kz+ky+2) = -trace(GT)*2*sigma_xi/(n*T);
SIG(kz+ky+2,kz+ky+1) = -trace(GT)*2*sigma_xi/(n*T);
SIG(kz+ky+3:kz+ky+2+J,kz+ky+1) = 0;  % These two last lines will change in terms of dimension
SIG(kz+ky+1,kz+ky+3:kz+ky+2+J) = 0;

% Elements corresponding to \sigma_xi\alpha, \sigma_xi\sigma_xi and \alpha\alpha
SIG(kz+ky+2,kz+ky+3:kz+ky+2+J) = 0;
SIG(kz+ky+3:kz+ky+2+J,kz+ky+2) = 0;
SIG(kz+ky+2,kz+ky+2) = 2*sigma_xi^2;

%sigma_e_sq = sigma_e*sigma_e;

for j1 = 1:J
    for j2 = 1:J
        SIG(kz+ky+2+j1,kz+ky+2+j2) = trace(sigma_e_sq*(sigma_e_inv==alpha(j1))*(sigma_e_inv==alpha(j2))) + ...
                          trace(sigma_e*(sigma_e_inv==alpha(j1))*sigma_e*(sigma_e_inv==alpha(j2))); 
    end
end

% Elements corresponding to \delta\beta_z
% These lines will change in terms of dimension
for k1 = 1:kz
    X_z_k1 = X_ZT(:,T*(k1-1)+1:k1*T);
    for l = 1:p
        SIG(k1,kz+ky+2+J+l) = -trace(M_Gy*kron(delta',eye(n))*X_z_k1*M_Fy*...
                               cZ(:,:,l)')/(n*T*sigma_xi^2);
        SIG(kz+ky+2+J+l,k1) = -trace(M_Gy*kron(delta',eye(n))*X_z_k1*M_Fy*...
                               cZ(:,:,l)')/(n*T*sigma_xi^2);   % cZ(:,:,l) = Z(:,:,l)-X_Z(:,:,l)*kron(beta_z(:,l),eye(T));
    end
end



% Elements corresponding to \delta\beta_y
% These lines will change in terms of dimension
for k1 = 1:ky
    X_y_k1 = X_Y(:,T*(k1-1)+1:k1*T);
    for l = 1:p
        SIG(kz+k1,kz+ky+2+J+l) = -trace(M_Gy*X_y_k1*M_Fy*cZ(:,:,l)')...
                                                        /(n*T*sigma_xi^2);
        SIG(kz+ky+2+J+l,kz+k1) = -trace(M_Gy*X_y_k1*M_Fy*cZ(:,:,l)')...
                                                        /(n*T*sigma_xi^2); % cZ(:,:,l) = Z(:,:,l)-X_Z(:,:,l)*kron(beta_z(:,l),eye(T));
    end
end

% Elements corresponding to \delta\lambda
% These lines will change in terms of dimension
for l = 1:p
    SIG(kz+ky+1,kz+ky+2+J+l) = -trace(M_Gy*X2*M_Fy*cZ(:,:,l)')...
                                                             /(n*T*sigma_xi^2);
    SIG(kz+ky+2+J+l,kz+ky+1) = -trace(M_Gy*X2*M_Fy*cZ(:,:,l)')...
                                                             /(n*T*sigma_xi^2);
end

% Elements corresponding to \delta\sigma_xi
% These lines will change in terms of dimension
for l = 1:p
    SIG(kz+ky+2,kz+ky+2+J+l) = 0;
    SIG(kz+ky+2+J+l,kz+ky+2) = 0;
end

% Elements corresponding to \delta\alpha
% These lines will change in terms of dimension
for l = 1:p
    for j = 1:J
        SIG(kz+ky+2+j,kz+ky+2+J+l) = 0;
        SIG(kz+ky+2+J+l,kz+ky+2+j) = 0;
    end
end
% Elements corresponding to \delta\delta
% These lines will change in terms of dimension
e = (1:p)';
for l_1 = 1:p
    for l_2 = 1:p
        SIG(kz+ky+2+J+l_1,kz+ky+2+J+l_2) = trace(M_Gy*cZ(:,:,l_1)*...
                       M_Fy*cZ(:,:,l_2)')/(n*T*sigma_xi^2) +...
                       trace((e==l_1)'*sigma_e_sq*(e==l_2))/sigma_xi^2; %cZ(:,:,l) = Z(:,:,l)-X_Z(:,:,l)*kron(beta_z(:,l),eye(T));
    end
end



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

bc = b_est - (C\phi)/sqrt(n*T);

% Theoretical standard deviation for the original estimator
%  obtained from diagonal elements of C^{-1}/nT.
CT = inv(C)/(n*T);
sd = sqrt(diag(CT));

% sd = zeros(6,1);
% sd(1) = sqrt(CT(1,1));
% sd(2) = sqrt(CT(2,2));
% sd(3) = sqrt(CT(3,3));
% sd(4) = sqrt(CT(4,4));
% sd(5) = sqrt(CT(5,5));
% sd(6) = sqrt(CT(6,6));

end
