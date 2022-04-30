% model_comparison_chapter5p6.m file
clear all;
sd = 221010;
rng(sd);
N = 10000;
xc = randn(N,1);
yc = randn(N,1);

Wmatrix(1).model = make_neighborsw(xc,yc,1);
for j=2:10
                    
        Wmatrix(j).model = make_neighborsw(xc,yc,j);

end

T = 10;
k = 2;
beta = [-1
        1];
theta = [0.5
         -0.5];    
bparms = [beta 
          theta];
x = randn(N*T,k);
sige = 0.001;
evec = randn(N*T,1)*sqrt(sige);

Wtrue = Wmatrix(5).model;
% SLX model x-matrix times beta
xbeta = [x  kron(eye(T),Wtrue)*x]*bparms;

% add fixed effects to the DGP
tts = (1:N)*(1/N);
SFE = kron(ones(T,1),tts');
ttt = (1:T)*(1/T);
TFE = kron(ttt',ones(N,1));
% slx model   
y = (xbeta + SFE + TFE + evec);

lmarginal_save = zeros(10,3);
rnames = strvcat('# of neighbors');
model = 3; % fixed effects for both regions and time periods
tic;
for iter = 1:10
    W = Wmatrix(iter).model;
    xmat = [x kron(eye(T),W)*x];

[ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,xmat,N,T,model);

info.lflag = 0; % exact log-determinant
% info.lflag = 1; % approximate log-determinant

result = lmarginal_static_panel(ywith,xwith,W,N,T,info); 

lmarginal_save(iter,:) = result.lmarginal';

rnames = strvcat(rnames,num2str(iter));

end
toc;

probs = model_probs(vec(lmarginal_save));

probs3 = reshape(probs,10,3);

in.rnames = rnames;
in.cnames = strvcat('SLX','SDM','SDEM');
in.width = 10000;
in.fmt = '%10.4f';
fprintf(1,'log-marginals for varying W-matrices AND models \n');
mprint(probs3,in);


