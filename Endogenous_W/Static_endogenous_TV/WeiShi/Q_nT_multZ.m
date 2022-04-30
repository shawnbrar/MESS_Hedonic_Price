% This function evaluates the likelihood function Q_nT(theta)
%   allow the z and y equations to have different number of factors
%
%   10/14/2015
% Modified: Cyrille D. 21/06/17

function f = Q_nT_multZ(theta,n,T,Rz,Ry,W,X_Y,Y,Z,X_Z,opt)

%global n T Rz Ry W X Y Z
% Determine 
%X_Z = X_z;
%X_Y = X;

kz = size(X_Z,2)/T;                               %the number of explanatory variables in Z equation
ky = size(X_Y,2)/T;                               %the number of explanatory variables in y equation
p  = size(Z,1)/n;                                %the number of dependent variables in the Z equation
J  = p + (p*(p-1))/2;                             %the number of distinct elements in Sigma_epsilon

beta_z      = theta(1:kz,1);
beta_y      = theta(kz+1:kz+ky,1);

if opt == 0
    lambda      = 0.9*(exp(2*theta(kz+ky+1))-1)/(exp(2*theta(kz+ky+1))+1);  % lambda < 0.9, > -0.9
    alpha_xi    = exp(theta(kz+ky+2));    % sigma_xi^{-1}
    alpha       = exp(theta(kz+ky+3:kz+ky+2+J));    % sigma_epsilon^{-1}
elseif opt == 1
    lambda      = theta(kz+ky+1);
    alpha_xi    = theta(kz+ky+2);
    alpha       = theta(kz+ky+3:kz+ky+2+J);
end
delta       = theta(kz+ky+2+J+1);


% alpha = (101)';
% p = 1;
sigma_eps = zeros(p,p);
lng = p-1:-1:0;
binf = 1;
for p_i = 1:p
    sigma_eps(p_i:p,p_i:p) = toeplitz(alpha(binf:binf+lng(p_i),1));
    binf = binf+lng(p_i) +1;
end

%sigma_eps_sq = 1/alpha^2; % replace by inv(Sigma_epsilon)
sigma_eps_sq = inv(sigma_eps.*sigma_eps);
sigma_xi_sq  = 1/alpha_xi^2;

ll = -(1/2)*log(det(sigma_eps_sq))-(1/2)*log(sigma_xi_sq);
for t = 1:T 
    S = eye(n)-lambda*W(:,:,t);
    ll = ll + log(det(S))/(n*T);
end

D = (kron(eye(n),inv(sigma_eps)))*(Z-X_Z*kron(beta_z,eye(T))); 
G = zeros(n,T);
cZ = Z-X_Z*kron(beta_z,eye(T));
cX_Y = X_Y*kron(beta_y,eye(T));

for t = 1:T
    S = eye(n)-lambda*W(:,:,t);
    G(:,t) = (1/sqrt(sigma_xi_sq))*(S*Y(:,t)-cX_Y(:,t)-kron(delta',eye(n))*(cZ(:,t)));
end

if sum(sum(isnan(D))) ~= 0 || sum(sum(isinf(D))) ~= 0 ||...
        sum(sum(isnan(G))) ~= 0 || sum(sum(isinf(G))) ~= 0 ||...
        isnan(ll) == 1 || isinf(ll) == 1
    disp('ll')
    disp(ll)
    disp('lambda')
    disp(lambda)
    disp('theta_3')
    disp(theta(3))
    disp('delta')
    disp(delta)
    disp('sigma_eps_sq')
    disp(sigma_eps_sq)
    disp('sigma_xi_sq')
    disp(sigma_xi_sq)
    f = 1/eps;
else
    [~,D1,~] = svd(D);
    D1 = diag(D1).*diag(D1);
    D1 = sort(D1);

    [~,G1,~] = svd(G);
    G1 = diag(G1).*diag(G1);
    G1 = sort(G1);

    ll = ll-(1/2)*sum(D1(1:end-Rz))/(n*T)-(1/2)*sum(G1(1:end-Ry))/(n*T);
    f = -ll;
end

end
