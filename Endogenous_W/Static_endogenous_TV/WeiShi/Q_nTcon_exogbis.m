% This function evaluates the likelihood function Q_nT(theta)
%   allow the z and y equations to have different number of factors
%
%   10/14/2015
% Created: Cyrille D. 20/10/17

function [fll,lly]  = Q_nTcon_exogbis(theta,n,T,Ry,W,x,y,info,biasopt)
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

ky = size(x,2)/T;  %the number of explanatory variables in y equation

beta_y      = theta(1:ky,1);

if strcmp(optimtyp,'Shi_con')==1
    lambda      = theta(ky+1);
    alpha_xi    = theta(ky+2);

    if boundtyp == 1
        sigma_xi_sq  = alpha_xi;
    elseif boundtyp == 2
        sigma_xi_sq  = 1/alpha_xi^2;
    end
elseif strcmp(optimtyp,'Shi_unc')==1
    if biasopt == 0
        lambda  = (exp(2*theta(ky+1))-1)/(exp(2*theta(ky+1))+1);  % lambda < 1, > -1
        alpha_xi= exp(theta(ky+2));    % sigma_xi^{-1}
   
        sigma_xi_sq  = 1/alpha_xi^2;
    elseif biasopt == 1
        lambda  = theta(ky+1);  
        alpha_xi= theta(ky+2);    
   
        sigma_xi_sq  = 1/alpha_xi^2;
    end
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
    f = 1/eps;
else

    [~,G1,~] = svd(G);
    G1 = diag(G1).*diag(G1);
    G1 = sort(G1);  
    
    ll = ll -(1/2)*sum(G1(1:end-Ry))/(n*T);
    fll = (-((1/2)*log(2*pi))+ll)*(n*T);
    lly = -(1/2)*sum(G1(1:end-Ry));
end

end
