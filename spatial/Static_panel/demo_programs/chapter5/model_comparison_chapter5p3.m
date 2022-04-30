% model_comparison_chapter5p2.m file
clear all;
sd = 221010;
rng(sd);

N = 400;
xc = randn(N,1);
yc = randn(N,1);

Wmatrix(1).model = make_neighborsw(xc,yc,1);

for j=2:10
                    
        Wmatrix(j).model = make_neighborsw(xc,yc,j);

end

T = 10;
k = 2;
beta = [1
        2];
theta = [1.5
         -0.5];    
bparms = [beta 
          theta];
x = randn(N*T,k);
rho = 0.4;
sige = 1;
evec = randn(N*T,1)*sqrt(sige);

Wtrue = Wmatrix(5).model;
xbeta = [x  kron(eye(T),Wtrue)*x]*bparms;

% add fixed effects to the DGP
tts = (1:N)*(1/N);
SFE = kron(ones(T,1),tts');
ttt = (1:T)*(1/T);
TFE = kron(ttt',ones(N,1));
    
y = (speye(N*T) - rho*kron(eye(T),Wtrue))\(xbeta + SFE + TFE + evec);

xmat = [x kron(eye(T),Wtrue)*x];

lmarginal_save = [];
rnames = strvcat('# of neighbors');

for iter = 1:10
    
    W = Wmatrix(iter).model;
    
model = 3; % fixed effects for both regions and time periods
[ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,xmat,N,T,model);

info.lflag = 0; % exact log-determinant
result = lmarginal_static_panel(ywith,xwith,W,N,T,info); 

lmarginal_save = [lmarginal_save
                  result.logm_sdm];
rnames = strvcat(rnames,num2str(iter));

end

probs = model_probs(lmarginal_save);

in.rnames = rnames;
in.cnames = strvcat('log-marginal','prob');
in.width = 10000;
in.fmt = '%10.4f';
fprintf(1,'log-marginals for varying W-matrices \n');
mprint([lmarginal_save probs],in);


