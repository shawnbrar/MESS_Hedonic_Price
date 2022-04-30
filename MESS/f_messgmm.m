function [f]=f_messgmm(par,y,x,M,Q,P,W)
[n k]=size(x);
b=par(2:end);
alp=par(1);

eps = expm(alp*M)*y-(x*b);
Weps=P*eps;
Z=[Weps Q];

for i = 1:n
m_t(i,:)=kron(eps(i,:),Z(i,:));
end
m=mean(m_t)';
f=m'*W*m;



