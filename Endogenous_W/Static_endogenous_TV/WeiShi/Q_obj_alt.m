function f = Q_obj_alt(theta,Y,W,Xy,Xz,Z,Ry,Rz)

%  In the paper, we need to compute the eigenvalues of DD' and not of D and the same for G (see lines 72 to 81). I
%  try with this new approach (DD') below
% Further, the eigenvalues are sorted in descending order and not in
% ascending order, since we need to consider the m largest 1. 
    % This function evaluates the likelihood function Q_nT(theta), see
    % Eq.6.
    % Ry(Rz) factors are concentrated out from y(z).
    N = (size(Y,1));
    T = (size(Y,2));
    Ky = (size(Xy,2)/T);
    Kz = (size(Xz,2)/T);

   
    beta_z  = theta(1:Kz,1);
    beta_y  = theta(Kz+1:Ky+Kz,1);
    lambda  = (exp(2*theta(Ky+Kz+1,1))-1)/(exp(2*theta(Ky+Kz+1,1))+1);  % lambda < 1, > -1
    alpha_xi= exp(theta(Ky+Kz+2,1));    % sigma_xi^{-1}
    alpha   = exp(theta(Ky+Kz+3,1));    % sigma_epsilon^{-1}
    delta   = theta(Ky+Kz+4,1);

    sigma_eps_sq = 1/alpha^2;
    sigma_xi_sq  = 1/alpha_xi^2;

    Xzsum = Xz*kron(beta_z,eye(T));
%     Xzsum = zeros(N,T);
%     for k = 1:Kz
%         Xzsum = Xzsum + beta_z(k)*squeeze(Xz(k,:,:));
%     end

    Xysum = Xy*kron(beta_y,eye(T));
%     Xysum = zeros(N,T);
%     for k = 1:Ky
%         Xysum = Xysum + beta_y(k)*squeeze(Xy(k,:,:));
%     end      
    
    ll = -(1/2)*log(sigma_eps_sq)-(1/2)*log(sigma_xi_sq);

    for t = 1:T 
        W1=squeeze(W(:,:,t));
        S = eye(N)-lambda*squeeze(W(:,:,t));
        ll = ll + log(det(S))/(N*T);
    end      
    
    D = (1/sqrt(sigma_eps_sq))*(Z-Xzsum);
    Ds=zeros(N,T);
    for t = 1 :T
        Ds(:,t) = (1/sqrt(sigma_eps_sq))*(Z(:,t)-Xzsum(:,t));
    end
    G = zeros(N,T);
    for t = 1:T
        S = eye(N)-lambda*squeeze(W(:,:,t));
        G(:,t) = (1/sqrt(sigma_xi_sq))*(S*Y(:,t)-Xysum(:,t)-(Z(:,t)-Xzsum(:,t))*delta);
    end    
    
if sum(sum(isnan(Ds))) ~= 0 || sum(sum(isinf(Ds))) ~= 0 || sum(sum(isnan(G))) ~= 0 || sum(sum(isinf(G))) ~= 0 || isnan(ll) == 1 || isinf(ll) == 1
    f = 1e15;
else
    [~,D1,~] = svd(D);
    D1 = diag(D1).*diag(D1);
    D1 = sort(D1);
    D2= sort(eig(Ds*Ds'),'descend');

    [~,G1,~] = svd(G);
    G1 = diag(G1).*diag(G1);
    G1 = sort(G1);
    G2=sort(eig(G*G'), 'descend');

    ll = ll-(1/2)*sum(D2(Rz+1:end))/(N*T)-(1/2)*sum(G2(Ry+1:end))/(N*T);
    f = -ll;
end    

end