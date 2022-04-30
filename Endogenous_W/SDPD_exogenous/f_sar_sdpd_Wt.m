function llike = f_sar_sdpd_Wt(parm,y,wy,x,W)

rho = 2/(1+exp(-parm))-1;

[n junk]=size(W);
[nt junk]=size(y);
[junk kx]=size(x);
t=nt/n; 

J=0;
for i=1:t
    St=speye(n)-rho*W(:,1+i*n:(i+1)*n);
    dSt=log(det(St));
    J=J+dSt;
end

Sy=y-rho*wy;
b = x\Sy;
e = Sy - x*b;
epe = e'*e;
%llike =(n/2)*log(epe/nt) -J/t;
llike =(n/2)*log(epe/nt) -J/t+(n/2)*(log(2*pi)+1);