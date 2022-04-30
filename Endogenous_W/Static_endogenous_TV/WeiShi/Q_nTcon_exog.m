% This function evaluates the likelihood function Q_nT(theta)
%   allow the z and y equations to have different number of factors
%
%   10/14/2015
% Created: Cyrille D. 20/10/17

function f = Q_nTcon_exog(theta,n,T,Ry,W,x,y,info)

    fieldsinf = fieldnames(info);
    nf = length(fieldsinf);
    for i=1:nf
        if strcmp(fieldsinf{i},'boundtyp')
            boundtyp = info.boundtyp;             
        end
    end
%global n T Rz Ry W X Y Z
% Determine 
%X_Z = X_z;
%X_Y = X;

ky = size(x,2)/T;  %the number of explanatory variables in y equation

beta_y      = theta(1:ky,1);
lambda      = theta(ky+1);
alpha_xi    = theta(ky+2);

if boundtyp == 1
    sigma_xi_sq  = alpha_xi;
elseif boundtyp == 2
    sigma_xi_sq  = 1/alpha_xi^2;
end

ll = -(1/2)*log(sigma_xi_sq);
for t = 1:T 
    S  = eye(n)-lambda*W(:,:,t);
    ll = ll + log(det(S))/(n*T);
end

G = zeros(n,T);
cX_Y = x*kron(beta_y,eye(T));

for t = 1:T
    S = eye(n)-lambda*W(:,:,t);
    G(:,t) = (1/sqrt(sigma_xi_sq))*(S*y(:,t)-cX_Y(:,t));
end

if sum(sum(isnan(G))) ~= 0 || sum(sum(isinf(G))) ~= 0 ||...
        isnan(ll) == 1 || isinf(ll) == 1
    disp('ll')
    disp(ll)
    disp('lambda')
    disp(lambda)
    disp('theta_3')
    disp(theta(3))
    disp('sigma_xi_sq')
    disp(sigma_xi_sq)
    f = 1/eps;
else

    [~,G1,~] = svd(G);
    G1 = diag(G1).*diag(G1);
    G1 = sort(G1);

    ll = ll -(1/2)*sum(G1(1:end-Ry))/(n*T);
    f = -ll;
end

end
