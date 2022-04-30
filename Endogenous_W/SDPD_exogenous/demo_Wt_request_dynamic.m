clear all;
%% Set sample size, seed and data
n = 100; 
T = 100;
nT=n*T;

for t=1:T 
coords=(rand(n,2)*10);
W{t}=make_neighborsw(coords(:,1), coords(:,2),5);
end
x=[ones(nT,1), randn(nT,2)];
bet=ones(3,1);
gam=0.4; %y_{t-1}
the=-0.5;% Wy_{t-1}
lam=0.4; %Wy_t

%%Generate the dependent variable where the initial values are set to 0.

Y=zeros(n*(T+1),1);
WY1=zeros(n*T,1);
Y1=zeros(n*T,1);

for t = 2:T
    A=eye(n)-lam*W{t};
Ai=eye(n)/A;
    stay=n*t+1;
    stopy=n*(t+1);
    stax=(t-1)*n + 1;
    stopx=n*t;
Y(stay:stopy,1)=Ai*(gam*Y(stax:stopx,1) + the*W{t-1}*Y(stax:stopx,1) +x(stax:stopx,:)*bet + randn(n,1));
WY1(stax:stopx,1)=W{t-1}*Y(stax:stopx,1);
Y1(stax:stopx,1)=Y(stax:stopx,1);
    
end

ymean=zeros(T,1);
for t = 1:T
     stax=(t-1)*n + 1;
    stopx=n*t;
    ymean(t,1) =mean(Y(stax:stopx,1));
end
plot(ymean)
xlabel('time') 
ylabel('y') 
title('Average of y for each period')
%% Estimation part    
% Remove the first 759 periods to get rid of the influence of the initial
% observations
% Y is of dimension n*(T+1) due to the initial period. 
T_start=50; % Date at which we start the estimation. 
x=x(end-(T_start*n)+1:end,:); %Create the X that will be used

% All variables are now vectors of 240 months *100 observations 

% Estimation of the SDPD model with time-constant matrices
Y_yu=Y(end-n*(T_start+1)+1:end); % Need to add the initial period because the code uses a (N*(T+1))vector for the dependent variable.
T1=length(Y_yu)/n;
info = struct('n',n,'t',T1,'rmin',-0.9,'rmax',0.9,'lflag',0,'tl',1,'stl',1); 

res_sdpd= sar_jihai(Y_yu,x(:,2:end),W{T-T_start+1},info); %Remove the constant because there are fixed effects and consider the initial weight matrix.
in1.rnames=strvcat('Variable','y_{t-1}','Wy_{t-1}','x1','x2','WY_t');
in1.cnames=strvcat(' Coeff.','t-stat');
disp('Estimation of the SDPD model with time-constant matrices')
out=[res_sdpd.beta res_sdpd.tstat2(1:end-2,1);
    res_sdpd.rho res_sdpd.tstat2(end-1,1)];
mprint(out,in1)


%%following is the dynamic spatial panel model
info = struct('n',n,'t',T1,'rmin',-0.9,'rmax',0.9,'lflag',0,'tl',1,'stl',1,'ted',1);
Wt=sparse(n,n*T_start+1);
for t=1:T_start+1
    beg=(t-1)*n+1;
    sto=(t)*n;
    Wt(:,beg:sto)=W{T-T_start-1+t}; % Attention Wt=[W1,W2, W3, ... WT]
end
result1 = sar_sdpd_Wt(Y_yu,x(:,2:end),Wt,info); %time varying Wt, See Lee and Yu (2012b)


% Use bias-corrected estimators 
out_2=[result1.theta1(1:end-2,1) result1.tstat1(1:end-2,1);
    result1.theta1(end-1,1) result1.tstat(end-1,1)];
disp('Estimation of the SDPD model with time-varying matrices')
mprint(out_2,in1)

true = [gam;
    the;
    bet(2:3,1);
    lam]

%y_t-1 = result1.beta(1,1)
%W_t-1y_{t-1}=result1.beta(2,1)
%X_t = result1.beta(3:end,1);
