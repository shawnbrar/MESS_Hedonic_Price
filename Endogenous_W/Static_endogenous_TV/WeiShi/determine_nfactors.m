function [Rhat] = determine_nfactors(est,W,y,x,z,x_z,equation,info)
    % Get informations needed for the estimation
    % Modified on 17/10/17 for constrained optimization (line 22-24)
    fieldsinf = fieldnames(info);
    nf = length(fieldsinf);
    for i=1:nf
        if strcmp(fieldsinf{i},'n')
            n = info.n;
        elseif strcmp(fieldsinf{i},'T')
            T = info.T;
        elseif strcmp(fieldsinf{i},'factors_choice_method')
            criteria = info.factors_choice_method;
        elseif strcmp(fieldsinf{i},'optimtyp')
            optimtyp = info.optimtyp; 
        elseif strcmp(fieldsinf{i},'boundtyp')
            boundtyp = info.boundtyp;             
        end
    end

    
    kz = size(x_z,2)/T;  %the number of explanatory variables in Z equation
    ky = size(x,2)/T;  %the number of explanatory variables in y equation

    % Get the estimated values of parameters
    beta_z      = est(1:kz,1);
    beta_y      = est(kz+1:kz+ky,1);  
    
    if strcmp(optimtyp,'Shi_con')==1
        lambda     =  est(kz+ky+1);
        alpha_xi   =  est(kz+ky+2);
        alpha      =  est(kz+ky+3);

        if boundtyp == 1
            sigma_eps_sq = alpha;
            sigma_xi_sq  = alpha_xi;
        elseif boundtyp == 2
            sigma_eps_sq = 1/alpha^2;
            sigma_xi_sq  = 1/alpha_xi^2;
        end
    elseif strcmp(optimtyp,'Shi_unc')==1
        lambda  = (exp(2*est(ky+kz+1))-1)/(exp(2*est(ky+kz+1))+1);  % lambda < 1, > -1
        alpha_xi= exp(est(ky+kz+2));    % sigma_xi^{-1}
        alpha   = exp(est(ky+kz+3));    % sigma_epsilon^{-1}
        
        sigma_eps_sq = 1/alpha^2;
        sigma_xi_sq  = 1/alpha_xi^2;
    end
    
    delta   = est(ky+kz+4);
    
    m      = min(n,T);
    neig_in_ratio = 21; % the number of eigenvalues to consider in the computation of ER and GR 
    
    if strcmp(equation,'endo_W') == 1
        D = (z-x_z*kron(beta_z,eye(T)))/sqrt(n*T*sigma_eps_sq);
        % list of (sorted) eigenvalues
        [~,D1,~] = svd(D);
        D1 = diag(D1).*diag(D1);  %squared to get the eigenvalues since svd gives the square root
        D1 = sort(D1,'descend');
        
        if strcmp(criteria,'Eigen_ratio')==1
            % Eigenvalue ratio criterion
            ER_z = zeros(neig_in_ratio,1);
            ER_z(1) = sum(D1)/log(m)/D1(1);  % Ratio for 0-th largest eigenvalue 
            for k = 2:neig_in_ratio
                ER_z(k) = D1(k)/D1(k+1);
            end
            
            Rhat    = find(ER_z==max(ER_z))-1;               % Why minus 1?
        elseif strcmp(criteria,'Growth_ratio')==1 
             % Growth ratio criterion
            GR_z = zeros(neig_in_ratio,1);
            for k = 1:neig_in_ratio
                if k == 1
                    GR_z(k) = log((sum(D1)+sum(D1)/log(m))/sum(D1(2:m)))/...
                              log(sum(D1(2:m))/sum(D1(3:m)));
                else
                    GR_z(k) = log((sum(D1(k:m)))/sum(D1(k+1:m)))/...
                              log(sum(D1(k+1:m))/sum(D1(k+2:m)));
                end
            end
%             [~,I]   = sort(GR_z, 'descend');
            Rhat   = find(GR_z==max(GR_z));
        end
    elseif strcmp(equation,'sar') == 1
        G = zeros(n,T);
        cz = z-x_z*kron(beta_z,eye(T));
        cx_y = x*kron(beta_y,eye(T));

        for t = 1:T
            %t =40;
            S = eye(n)-lambda*W(:,:,t);
            G(:,t) = S*y(:,t)-cx_y(:,t)-(cz(:,t))*delta;
        end
        G = G/sqrt(n*T*sigma_xi_sq);
        
        [~,G1,~] = svd(G);
        G1 = diag(G1).*diag(G1);
        G1 = sort(G1,'descend');
        
        if strcmp(criteria,'Eigen_ratio')==1
            % Eigenvalue ratio criterion
            ER_y = zeros(neig_in_ratio,1);
            ER_y(1) = sum(G1)/log(m)/G1(1);
            for k = 2:neig_in_ratio
                ER_y(k) = G1(k)/G1(k+1);
            end
           
            Rhat   = find(ER_y==max(ER_y))-1;
           
        elseif strcmp(criteria,'Growth_ratio')==1 
            % Growth ratio criterion
            GR_y = zeros(neig_in_ratio,1);
            for k = 1:neig_in_ratio
                if k == 1
                    GR_y(k) = log((sum(G1)+sum(G1)/log(m))/sum(G1(2:m)))/...
                              log(sum(G1(2:m))/sum(G1(3:m)));
                else
                    GR_y(k) = log((sum(G1(k:m)))/sum(G1(k+1:m)))/...
                              log(sum(G1(k+1:m))/sum(G1(k+2:m)));
                end
            end
           Rhat   = find(GR_y==max(GR_y));
        end
    end
    if Rhat >=T
        Rhat=T-1;
    end
end