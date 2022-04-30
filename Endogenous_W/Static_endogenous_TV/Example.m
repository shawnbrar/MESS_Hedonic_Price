
% Create a static panel data model with time-varying endogenous W and
% factors
clear

 n = 50;
 T=20;
 Sigma=[1 0.5;
     0.5 2];
 rh=0.4;
 x1=randn(n*T,2);
 x2=randn(n*T,2);
 bety=[1;1];
 betz=[0.5;0.5];
 In=eye(n);
 err=mvnrnd([0,0],Sigma,n*T);
 % Create factors
Rz0 = 1;      % true number of factors
Ry0 = 4;      % true number of factors in the y equation

Gamma_y = rand(n,Ry0)*4-ones(n,Ry0)*2;
Gamma_z = rand(n,Rz0)*4-ones(n,Rz0)*2;
F_y     = rand(T,Ry0)*4-ones(T,Ry0)*2;
F_z     = rand(T,Rz0)*4-ones(T,Rz0)*2;
y=zeros(n*T,1);
 for t = 1:T
     
     b=(t-1)*n+1;
     e=n*t;
     Z(:,t)=x2(b:e,:)*betz+Gamma_z*F_z(t,1)+err(b:e,2);
     for i = 1:n
            for j = 1:n
                W(i,j,t) =1/(1+abs(Z(i,t)-Z(j,t)));
            end
     end
     W(:,:,t)=normw(W(:,:,t));
     Ai=In/(In-rh*W(:,:,t));
     
     Fac_y=zeros(n,1);
     for f=1:Ry0
         Fac_y=Fac_y + Gamma_y(:,f)*F_y(t,f);
     end
     y(b:e,1)=Ai*(x1(b:e,:)*bety +Fac_y +err(b:e,1));
 end
 Z_fin=reshape(Z,n*t,1);
 Rys=10;
 Rzs=10;
 optimtyp ='Shi_unc';
 info = struct('n',n,'T',t,'Ry',Rys,...
                        'Rz',Rzs, 'optimtyp',optimtyp,...
                        'TolFun',1e-10,...
                        'TolX',1e-10,'factors_choice',1,...
                  'factors_choice_method',{'Growth_ratio'},'cbiais',1,...
                  'print_res',1);
  [outputs,allres] = launch_estimation_shi(y,x1,Z_fin, [x1,x2],W,info);
 kz = size([x1,x2],2);
  outputs.allestimate(1:kz)  %Estimates for the regressors entering
  %   equation for z
   ky = size(x1,2) 
outputs.allestimate(1+kz:kz+ky) %Estimates for the regressod entering the
    % main equation

outputs.allestimate(kz+ky+1,1) %Estimates for Wy
outputs.allestimate(end,1) %Estimates for delta (assessing the role of
% endogeneity

outputs.allnfactors% Number of factors for y and z (y first and z second)




 
