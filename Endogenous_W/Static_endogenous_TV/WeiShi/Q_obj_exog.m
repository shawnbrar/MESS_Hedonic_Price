function f = Q_obj_exog(theta,Y,W,Xy,Ry)
    % This function evaluates the likelihood function Q_nT(theta), see
    % Eq.6.
    % Ry(Rz) factors are concentrated out from y(z).
    N = size(Y,1);
    T = size(Y,2);
    Ky = size(Xy,2)/T;
   
    
    beta_y  = theta(1:Ky);
    lambda  = (exp(2*theta(Ky+1))-1)/(exp(2*theta(Ky+1))+1);  % lambda < 1, > -1
    alpha_xi= exp(theta(Ky+2));    % sigma_xi^{-1}

    sigma_xi_sq  = 1/alpha_xi^2;

    Xysum = Xy*kron(beta_y,eye(T));
%     Xysum = zeros(N,T);
%     for k = 1:Ky
%         Xysum = Xysum + beta_y(k)*squeeze(Xy(k,:,:));
%     end      
    
    ll = -(1/2)*log(sigma_xi_sq);
    for t = 1:T 
        S = eye(N)-lambda*squeeze(W(:,:,t));
        ll = ll + log(det(S))/(N*T);
    end      
    
    G = zeros(N,T);
    for t = 1:T
        S = eye(N)-lambda*squeeze(W(:,:,t));
        G(:,t) = (1/sqrt(sigma_xi_sq))*(S*Y(:,t)-Xysum(:,t));
    end    
    
if  sum(sum(isnan(G))) ~= 0 || sum(sum(isinf(G))) ~= 0 ...
        || isnan(ll) == 1 || isinf(ll) == 1
    f = 4.503e+15;
else
  
    [~,G1,~] = svd(G);
    G1 = diag(G1).*diag(G1);
    G1 = sort(G1);

    ll = ll - (1/2)*sum(G1(1:end-Ry))/(N*T);
    f = -ll;
end    

end