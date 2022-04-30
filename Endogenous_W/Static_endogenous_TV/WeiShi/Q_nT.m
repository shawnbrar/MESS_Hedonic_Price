% This function evaluates the likelihood function Q_nT(theta)
%   allow the z and y equations to have different number of factors
%
%   10/14/2015
% Modified: Cyrille D. 21/06/17

function f = Q_nT(theta,n,T,Rz,Ry,W,x,y,z,x_z,opts)

%global n T Rz Ry W X Y Z
% Determine 
%X_Z = X_z;
%X_Y = X;

kz = size(x_z,2)/T;  %the number of explanatory variables in Z equation
ky = size(x,2)/T;  %the number of explanatory variables in y equation

beta_z      = theta(1:kz,1);
beta_y      = theta(kz+1:kz+ky,1);

if opts == 0
    lambda      = 0.9*(exp(2*theta(kz+ky+1))-1)/(exp(2*theta(kz+ky+1))+1);  % lambda < 0.9, > -0.9
    alpha_xi    = exp(theta(kz+ky+2));    % sigma_xi^{-1}
    alpha       = exp(theta(kz+ky+3));    % sigma_epsilon^{-1}
elseif opts == 1
    lambda      = theta(kz+ky+1);
    alpha_xi    = theta(kz+ky+2);
    alpha       = theta(kz+ky+3);
end
delta       = theta(kz+ky+4);

sigma_eps_sq = 1/alpha^2;
sigma_xi_sq  = 1/alpha_xi^2;

ll = -(1/2)*log(sigma_eps_sq)-(1/2)*log(sigma_xi_sq);
for t = 1:T 
    S = eye(n)-lambda*W(:,:,t);
    ll = ll + log(det(S))/(n*T);
end

D = (1/sqrt(sigma_eps_sq))*(z-x_z*kron(beta_z,eye(T)));       % Include a Kronecker product to allow z to be of dimension p
G = zeros(n,T);
cZ = z-x_z*kron(beta_z,eye(T));
cX_Y = x*kron(beta_y,eye(T));

for t = 1:T
    S = eye(n)-lambda*W(:,:,t);
    G(:,t) = (1/sqrt(sigma_xi_sq))*(S*y(:,t)-cX_Y(:,t)-(cZ(:,t))*delta);
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
