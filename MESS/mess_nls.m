function results=mess_nls(y,x,W)
% Purpose: estimation  of the MESS models in simulations
amin=-5;
amax=5;
parm=0; %Initial value for alpha
% if nargin==4
%      if ~isstruct(info)
%  error('sarar: must supply the options as a structure variable');
%  end;
% info.MaxIter = 500;
%  fields = fieldnames(info);
%  nf = length(fields);
%  for i=1:nf
%     if strcmp(fields{i},'parm')
%        parm = info.parm;
% 	elseif strcmp(fields{i},'convg')
%        options.TolFun = info.convg;
%     elseif strcmp(fields{i},'maxit')
%         options.MaxIter  = info.maxit;
%     elseif strcmp(fields{i},'Nhes')
%         Nhes  = info.Nhes;
%     elseif strcmp(fields{i},'order')
%         order = info.order;  results.order = order;
%     elseif strcmp(fields{i},'iter')
%     liter = info.iter; results.liter = liter;
%     end;
%  end;
% else nargin == 3 % use default options
% options = optimset('fminsearch');
% end

[n, k]=size(x);
% options=optimset('OutputFcn',@outfun);
options=optimset('fminsearch');
options= optimset(options, 'TolX',10^(-6));
f_mess_oc = @(r) f_mess(r,y,x,W);
[alp fval exit output] = fminsearch(f_mess_oc,parm,options);
results.exit=exit;
results.alpha=alp;
ay=expm(alp*W)*y;

b=(x'*x)\(x'*ay);
%%%%%%% Use b=(x'*x)\(x'*ay) to save a little time, as x'*ay is a vector. 
results.beta=b;
res=ay-x*b;
results.resid=res;
sige=res'*res/n;
results.sigma2=sige;
par=[results.beta
    results.alpha
    results.sigma2];
% results.lik=f_mess2(par,y,x,W);

% % Computation of the information matrix (for Normal and homoskedastic
% % errors)

Wxb=W*(x*b);
%%%%%%% Use Wxb=W*(x*b)
xpx = zeros(k+2);
xpx(1:k,1:k)=1/sige*x'*x;
xpx(1:k,k+1)=-1/sige*x'*Wxb;
xpx(k+1,1:k)= xpx(1:k,k+1)';
xpx(k+2,k+2)=n/(2*sige^2);
%%%%%%%%% It should be xpx(k+2,k+2)=n/(2*sige^2), so missing 2.
xpx(k+1,k+1)=sum(sum(W.*(W'+W)))+1/sige*Wxb'*Wxb;
%%%%%% Using sum(sum(W.*(W'+W))) is faster than trace(W'*(W+W')), since
%%%%%% sum(sum(A'.*B)) computes trace(A*B).

xpxi=eye(k+2)/xpx;
temp=diag(xpxi(1:k+1,1:k+1));
results.cov = xpxi(1:end-1,1:end-1);
results.tstat=[results.beta; results.alpha]./sqrt(temp);
% 
%  dhessn = hessian('f_mess2',par,y,x,W);
% hessi = invpd(-dhessn);

% % Computation of the sandwich form for QML with heteroskedastic errors
Im=zeros(k+1);
% Creation of the heteroskedastic errors matrix
Wxb=W*(x*b);
s=res.*res;
Sigm=diag(s);
results.Sigm = Sigm;
WpW=W'+W;
WpWpW=W'.*WpW;
tr=sum(WpWpW,2);
TR=sum(s.*tr,1);
%%%%%%%% I am not sure whether line 84 is right,  but another way to
%%%%%%%% compute is TR=sum(sum((sparse(Sigm)*W').*(W'+W))), I use a
%%%%%%%% sparse(Sigm) so that Sigma*W' can be computed quickly. Note
%%%%%%%% trace(A*B)=trace(A'*B) if B is symmetric. It works but it's slower
%%%%%%%% than what is used
% Information matrix
Im(1:k,1:k)=x'*x;
Im(1:k,k+1)=-x'*Wxb;
Im(k+1,1:k)=Im(1:k,k+1)';
Im(k+1,k+1)=TR+Wxb'*Wxb;
% Inverse of the information matrix
Ii=eye(k+1)/(Im);
% Outer product of the score

TR2=trace(Sigm*W'*Sigm*(WpW));
Sc=zeros(k+1);
Sc(1:k,1:k)=x'*Sigm*x;
Sc(1:k,k+1)=-x'*Sigm*Wxb;
%%%%%%%%%%% Should be Sc(1:k,k+1)=-x'*Sigm*Wxb, missing "-".
Sc(k+1,1:k)=Sc(1:k,k+1)';
Sc(k+1,k+1)=TR2+Wxb'*Sigm*Wxb;
xpxi=Ii*Sc*Ii;
temp=diag(xpxi);
results.covr = xpxi;
results.tstatr=[results.beta;results.alpha]./sqrt(temp);
parh=par(1:end-1); %remove sigma from the parameter vector since it does not exist in heteroskedastik case
results.likr=f_mess2he(parh, y,x,W,Sigm);

function llik=f_mess2he(parm,y,x,W,Sigm)
% Purpose = concentrated function to minimize over in order to get an
% estimated value of alpha

% y is the vector of the dependent variable (N x 1)
% x is the  matrix of independent variables (N x k)
% W is the weight matrix
[n,k]=size(x);
b=parm(1:k);
alp=parm(k+1);
res=expm(alp*W)*y-x*b;
Si=inv(Sigm);
epe=res'*Si*res;
ds=diag(Sigm);
lds=log(ds);
ilds = ones(1,n)*lds;
llik=-n/2*(log(2*pi))-ilds/2-1/(2)*epe;
