% This function evaluates the likelihood function Q_nT(theta)
%   allow the z and y equations to have different number of factors
%
%   10/14/2015
% Modified: Cyrille D. 21/06/17

%  In the paper, we need to compute the eigenvalues of DD' and not of D and the same for G (see lines 72 to 81). I
%  try with this new approach (DD') below
% Further, the eigenvalues are sorted in descending order and not in
% ascending order, since we need to consider the m largest 1. 

function [f] = Q_nTcon_alt(theta,n,T,Rz,Ry,W,x,y,z,x_z,info)
    fieldsinf = fieldnames(info);
    nf = length(fieldsinf);
    for i=1:nf
        if strcmp(fieldsinf{i},'boundtyp')
            boundtyp = info.boundtyp;             
        end
    end
% Determine 
%X_Z = X_z;
%X_Y = X;
%theta = initial_values;

kz = size(x_z,2)/T;  %the number of explanatory variables in Z equation
ky = size(x,2)/T;  %the number of explanatory variables in y equation

beta_z      = theta(1:kz,1);
beta_y      = theta(kz+1:kz+ky,1);
lambda      = theta(kz+ky+1,1);
alpha_xi    = theta(kz+ky+2,1); %sigma_xi^2
alpha       = theta(kz+ky+3,1); %sigma_eps^2
delta       = theta(kz+ky+4,1);

if boundtyp == 1
    sigma_eps_sq = alpha;
    sigma_xi_sq  = alpha_xi;
elseif boundtyp == 2
    sigma_eps_sq = 1/alpha^2;
    sigma_xi_sq  = 1/alpha_xi^2;
end

cX_Z = x_z*kron(beta_z,eye(T));
cX_Y = x*kron(beta_y,eye(T));

ll = -(1/2)*log(sigma_eps_sq)-(1/2)*log(sigma_xi_sq);
for t = 1:T 
    S = eye(n)-lambda*squeeze(W(:,:,t));
    ll = ll + log(det(S))/(n*T);
end


D = (1/sqrt(sigma_eps_sq))*(z-cX_Z);% Include a Kronecker product to allow z to be of dimension p
% 
Ds=zeros(n,T);
for t=1:t
  Ds(:,t)= (1/sqrt(sigma_eps_sq))*(z(:,t)-cX_Z(:,t));
end

DDp=Ds*Ds';


G = zeros(n,T);
for t = 1:T
    S = eye(n)-lambda*W(:,:,t);
    G(:,t) = (1/sqrt(sigma_xi_sq))*(S*y(:,t)-cX_Y(:,t)-(z(:,t) - cX_Z(:,t))*delta);
end

if sum(sum(isnan(Ds))) ~= 0 || sum(sum(isinf(Ds))) ~= 0 ||...
        sum(sum(isnan(G))) ~= 0 || sum(sum(isinf(G))) ~= 0 ||...
        isnan(ll) == 1 || isinf(ll) == 1
 
    f = 1/eps;
else
    [~,D1,~] = svd(D);
    D1 = diag(D1).*diag(D1);
%     D1 = sort(D1);
    D2 = sort(eig(DDp),'descend');
    

    [~,G1,~] = svd(G);
    G1 = diag(G1).*diag(G1);
    G1 = sort(G1);
    G2= sort(eig(G*G'),'descend');

    
    ll2=ll-(1/2)*sum(D2(Rz+1:end))/(n*T)-(1/2)*sum(G2(Ry+1:end))/(n*T);
    f=-ll2;
    
end

end
