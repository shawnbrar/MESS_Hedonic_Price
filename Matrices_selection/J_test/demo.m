n = 300;
coords=rand(n,2);
Wt=make_neighborsw(coords(:,1),coords(:,2),5);
x=[ones(n,1), randn(n,2)];
bet=[1;2;3];
xb=x*bet;
e=randraw('t',[10 0 1],n,1).*rand(n,1);
%Construction of the reduced form of the SAR
In=eye(n);
rh=0.5;
A=In-rh*Wt;
Ai=In/A;
y=Ai*(xb + e);
Wx=Wt*x(:,2:end);
W2x=Wt*Wx;
Q=[x, Wx W2x];
P=Wt;
res_1=gmmestimation_sar_he(y,x,Wt,Q,P);
coords1=rand(n,2);
W2=dinvarc(coords1(:,1),coords1(:,2),1);
W2=normw(W2);
P2=W2;
res_2=gmmestimation_sar_he(y,x,W2,Q,P2);
coords2=rand(n,2);
W3=make_neighborsw(coords2(:,1),coords2(:,2),5);
P3=W3;
res_3=gmmestimation_sar_he(y,x,W3,Q,P3);

B=59; % Number of bootstrap.
[stat11, stat21,stat31,signi1]=Jp1_hetero_3(y,Wt,W2,W3,x,res_1,res_2,res_3,B)
[stat12, stat22,stat32,signi2]=Jp2_hetero_3(y,Wt,W2,W3,x,res_1,res_2,res_3,B)


