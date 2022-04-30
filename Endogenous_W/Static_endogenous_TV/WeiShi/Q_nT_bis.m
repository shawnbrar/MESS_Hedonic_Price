% This function evaluates the likelihood function Q_nT(theta)
%   allow the z and y equations to have different number of factors
%
%   10/14/2015
% Modified: Cyrille D. 21/06/17

function [fll,lly,llz] = Q_nT_bis_alt(theta,n,T,Rz,Ry,W,X_Y,Y,Z,X_Z,info,biasopt)
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
% Determine 
%X_Z = X_z;
%X_Y = X;

kz = size(X_Z,2)/T;  %the number of explanatory variables in Z equation
ky = size(X_Y,2)/T;  %the number of explanatory variables in y equation
p  = size(Z,3);                    %the number of dependent variables in the Z equation

beta_z      = theta(1:kz,1);
beta_y      = theta(kz+1:kz+ky,1);

if strcmp(optimtyp,'Shi_con')==1
    lambda      = theta(kz+ky+1);
    alpha_xi    = theta(kz+ky+2);
    alpha       = theta(kz+ky+3);

    if boundtyp == 1
        sigma_eps_sq = alpha;
        sigma_xi_sq  = alpha_xi;
    elseif boundtyp == 2
        sigma_eps_sq = 1/alpha^2;
        sigma_xi_sq  = 1/alpha_xi^2;
    end
elseif strcmp(optimtyp,'Shi_unc')==1
    if biasopt == 0
        lambda  = (exp(2*theta(ky+kz+1))-1)/(exp(2*theta(ky+kz+1))+1);  % lambda < 1, > -1
        alpha_xi= exp(theta(ky+kz+2));    % sigma_xi^{-1}
        alpha   = exp(theta(ky+kz+3));    % sigma_epsilon^{-1}

        sigma_eps_sq = 1/alpha^2;
        sigma_xi_sq  = 1/alpha_xi^2;
    elseif biasopt == 1
        lambda  = theta(ky+kz+1);  
        alpha_xi= theta(ky+kz+2);    
        alpha   = theta(ky+kz+3);   

        sigma_eps_sq = 1/alpha^2;
        sigma_xi_sq  = 1/alpha_xi^2;
    end
end

delta       = theta(kz+ky+4);


ll = -(1/2)*log(sigma_eps_sq)-(1/2)*log(sigma_xi_sq);
for t = 1:T 
    S = eye(n)-lambda*W(:,:,t);
    ll = ll + log(det(S))/(n*T);
end

D = (1/sqrt(sigma_eps_sq))*(Z-X_Z*kron(beta_z,eye(T)));       % Include a Kronecker product to allow z to be of dimension p
G = zeros(n,T);
cZ = Z-X_Z*kron(beta_z,eye(T));
cX_Y = X_Y*kron(beta_y,eye(T));

for t = 1:T
    S = eye(n)-lambda*W(:,:,t);
    G(:,t) = (1/sqrt(sigma_xi_sq))*(S*Y(:,t)-cX_Y(:,t)-(cZ(:,t))*delta);
end

if sum(sum(isnan(D))) ~= 0 || sum(sum(isinf(D))) ~= 0 ||...
        sum(sum(isnan(G))) ~= 0 || sum(sum(isinf(G))) ~= 0 ||...
        isnan(ll) == 1 || isinf(ll) == 1

    f = 1/eps;
else
    
    
    [~,D1,~] = svd(D);
    D1 = diag(D1).*diag(D1);
    D1 = sort(D1,'descend');
%     sD = size(D)
%     sD1 = size(D1)

    [~,G1,~] = svd(G);
    G1 = diag(G1).*diag(G1);
    G1 = sort(G1);
%     sG = size(G)
%     sG1 = size(G1)
    
    
    ll = ll-(1/2)*sum(D1(1:end-Rz))/(n*T)-(1/2)*sum(G1(1:end-Ry))/(n*T);
    fll = (-(((p+1)/2)*log(2*pi))+ll)*(n*T);
    llz = -(1/2)*sum(D1(1:end-Rz));
    lly = -(1/2)*sum(G1(1:end-Ry));
%     
%     zD1 = D1(1:end-Rz)
%     yG1 = G1(1:end-Ry)
end
end
