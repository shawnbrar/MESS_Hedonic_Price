% This function evaluates the likelihood function Q_nT(theta)
%   ignoring the endogeneity of the weights matrix
%   The factors are considered.
%
%   10/15/2015

function f = Q_nT_SAR(theta)

global n T Ry W X Y

beta_y      = theta(1);
lambda      = 0.9*(exp(2*theta(2))-1)/(exp(2*theta(2))+1);  % lambda < 0.9, > -0.9
alpha_xi    = exp(theta(3));                                % sigma_xi^{-1}

sigma_xi_sq  = 1/alpha_xi^2;

ll = -(1/2)*log(sigma_xi_sq);
for t = 1:T 
    S = eye(n)-lambda*W(:,:,t);
    ll = ll + log(det(S))/(n*T);
end

G = zeros(n,T);
for t = 1:T
    S = eye(n)-lambda*W(:,:,t);
    G(:,t) = (1/sqrt(sigma_xi_sq))*(S*Y(:,t)-X(:,t)*beta_y);
end

if sum(sum(isnan(G))) ~= 0 || sum(sum(isinf(G))) ~= 0 || isnan(ll) == 1 || isinf(ll) == 1
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

    ll = ll-(1/2)*sum(G1(1:end-Ry))/(n*T);
    f = -ll;
end

end
