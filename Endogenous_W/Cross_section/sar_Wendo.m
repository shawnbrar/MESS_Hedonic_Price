function results = sar_Wendo(y,x1,W,Zn,x2,info)
% PURPOSE: computes spatial autoregressive model estimates
%           y = p*W*y + X*b + e, using sparse matrix algorithms
% ---------------------------------------------------
%  USAGE: results = sar(y,x,W,info)
%  where:  y = dependent variable vector
%          x = explanatory variables matrix, (with intercept term in first
%              column if used)
%          W = standardized weight matrix 
%       info = an (optional) structure variable with input options:
%       info.rmin  = (optional) minimum value of rho to use in search (default = -1) 
%       info.rmax  = (optional) maximum value of rho to use in search (default = +1)    
%       info.eig   = 0 for default rmin = -1,rmax = +1, 1 for eigenvalue calculation of these
%       info.convg = (optional) convergence criterion (default = 1e-8)
%       info.maxit = (optional) maximum # of iterations (default = 500)
%       info.lflag = 0 for full lndet computation (default = 1, fastest)
%                  = 1 for MC lndet approximation (fast for very large problems)
%                  = 2 for Spline lndet approximation (medium speed)
%       info.order = order to use with info.lflag = 1 option (default = 50)
%       info.iter  = iterations to use with info.lflag = 1 option (default = 30)  
%       info.lndet = a matrix returned by sar, sar_g, sarp_g, etc.
%                    containing log-determinant information to save time
%       info.ndraw = 1,000 by default
% ---------------------------------------------------
%  RETURNS: a structure
%         results.meth  = 'sar'
%         results.beta  = bhat (nvar x 1) vector
%         results.rho   = rho
%         results.tstat = asymp t-stat (last entry is rho)
%         results.bstd  = std of betas (nvar x 1) vector
%         results.pstd  = std of rho
%         results.total    = a 3-d matrix (ndraw,p,ntrs) total x-impacts
%         results.direct   = a 3-d matrix (ndraw,p,ntrs) direct x-impacts
%         results.indirect = a 3-d matrix (ndraw,p,ntrs) indirect x-impacts
%                            ndraw = 1,000 by default, ntrs = 101 default
%                            p = nvar-1 if there is a constant term which we skip
%         results.yhat  = yhat         (nobs x 1) vector
%         results.resid = residuals    (nobs x 1) vector
%         results.sige  = sige = (y-p*W*y-x*b)'*(y-p*W*y-x*b)/n
%         results.rsqr  = rsquared
%         results.rbar  = rbar-squared
%         results.lik   = log likelihood
%         results.nobs  = # of observations
%         results.nvar  = # of explanatory variables in x 
%         results.p     = # of explanatory variables excluding the constant term
%         results.cflag = 0 for no intercept term, 1 for intercept term
%         results.y     = y data vector
%         results.iter  = # of iterations taken
%         results.rmax  = 1/max eigenvalue of W (or rmax if input)
%         results.rmin  = 1/min eigenvalue of W (or rmin if input)
%         results.lflag = lflag from input
%         results.liter = info.iter option from input
%         results.order = info.order option from input
%         results.limit = matrix of [rho lower95,logdet approx, upper95] intervals
%                         for the case of lflag = 1
%         results.time1 = time for log determinant calcluation
%         results.time2 = time for eigenvalue calculation
%         results.time3 = time for hessian or information matrix calculation
%         results.time4 = time for optimization
%         results.time5 = time for effects estimates calculation
%         results.time  = total time taken      
%         results.lndet = a matrix containing log-determinant information
%                          (for use in later function calls to save time)
%         results.ndraw = # of draws used for x-impact estimates (1,000 by default)
% --------------------------------------------------
%  NOTES: if you use lflag = 1 or 2, info.rmin will be set = -1 
%                                    info.rmax will be set = 1
%         For n < 500 you should use lflag = 0 to get exact results                                    
% --------------------------------------------------  
%  SEE ALSO: sar_d, sar_d2 demos
% ---------------------------------------------------
% REFERENCES: LeSage and Pace (2009) Chapter 4 on maximum likelihood estimation 
%             of spatial regression models.
% For lndet information see: Chapter 4
% For interpretation of direct, indirect and total x-impacts see: Chapter 2
% ---------------------------------------------------

% written by:
% James P. LeSage, last modified 3/2010
% Dept of Finance & Economics
% Texas State University-San Marcos
% 601 University Drive
% San Marcos, TX 78666
% jlesage@spatial-econometrics.com


time1 = 0; 
time2 = 0;
time3 = 0;
time4 = 0;
time5 = 0;

timet = clock; % start the clock for overall timing

% check size of user inputs for comformability
[n, nvar] = size(x1); 
[n1, n2] = size(W);

% if we have no options, invoke defaults
if nargin == 3
    info.lflag = 1;
end;


% parse input options
[rmin,rmax,convg,maxit,detval,ldetflag,eflag,order,miter,options,ndraw] = sar_parse(info);

results.ndraw = ndraw;
    
% check if the user handled the intercept term okay
    n = length(y);
    if sum(x1(:,1)) ~= n
    tst = sum(x1); % we may have no intercept term
    ind = find(tst == n); % we do have an intercept term
     if length(ind) > 0
     error('sar: intercept term must be in first column of the x-matrix');
     elseif length(ind) == 0 % case of no intercept term
     cflag = 0;
     p = size(x1,2);
     end;
    elseif sum(x1(:,1)) == n % we have an intercept in the right place
     cflag = 1;
     p = size(x1,2)-1;
    end;
     
    results.cflag = cflag;
    results.p = p;
    
    if n1 ~= n2
    error('sar: wrong size weight matrix W');
    elseif n1 ~= n
    error('sar: wrong size weight matrix W');
    end;
    [nchk junk] = size(y);
    if nchk ~= n
    error('sar: wrong size y vector input');
    end;
    
 

% compute eigenvalues or limits
[rmin,rmax,time2] = sar_eigs(eflag,W,rmin,rmax,n);
results.time2 = time2;

% do log-det calculations
[detval,time1] = sar_lndet(ldetflag,W,rmin,rmax,detval,order,miter);
results.time1 = time1;




%%

t0 = clock;



% step 2) maximize concentrated likelihood function for thetahat;
options=optimset('fmincon');
options = optimset(options,'Display', 'off');
options = optimset(options,'FunValCheck', 'on');
%options = optimset(options,'Diagnostics', 'off');
options.Algorithm = 'interior-point';
options = optimset(options, 'TolFun',10^(-05));
par0=[0;zeros(size(x1,2),1);zeros(size(x2,2),1);0.5;0.5;0.5];

%%%%%
%Thea=[alpha,beta,gamma,sig2_ep,sig2_v,coeffcorr]
lb=[rmin+0.001;-ones(size(x1,2),1)*Inf;-ones(size(x2,2),1)*Inf;0.001;0.001;-1]; % -1 is the minimum value for lambda and 0 for theta because the function is symmetric for theta. We can thus maximize over [0, infinity]  
ub=[rmax-0.001;ones(size(x1,2),1)*Inf;ones(size(x2,2),1)*Inf;Inf;Inf;1];



f_sar_Wendo_oc=@(r) f_sar_Wendo(r,detval,n,Zn,x2,y,W,x1);
[ptheta,liktmp,exitflag,output] = fmincon(f_sar_Wendo_oc,par0,[],[],[],[],lb,ub,[],options);
%[ptheta,liktmp,exitflag,output] = fminunc('f_sar_Wendo',par0,options,detval,n,Zn,x2,y,W,x1); 


    
time4 = etime(clock,t0);

%if exitflag == 0
%fprintf(1,'\n sar: convergence not obtained in %4d iterations \n',output.iterations);
%end;
%results.iter = output.iterations;


%%


prho=ptheta(1,1);
pbeta=ptheta(2:size(x1,2)+1,1);
pGamma=ptheta(size(x1,2)+2:size(x1,2)+1+size(x2,2),1); %hypoth�se que p2=1
psig2_ep=ptheta(size(x1,2)+2+size(x2,2),1); %hypoth�se que p2=1
psig2_v=ptheta(size(x1,2)+2+size(x2,2)+1,1);
pcoeffcorr=ptheta(size(x1,2)+2+size(x2,2)+2,1); %hypoth�se que p2=1
%calcul suppl�mentaire en partant des estimations
psig_vep=pcoeffcorr*sqrt(psig2_v)*sqrt(psig2_ep);
pdelta=inv(psig2_ep)*psig_vep;
psig2_eta=psig2_v-psig_vep'*inv(psig2_ep)*psig_vep;



% step 3) find b,sige maximum likelihood estimates
results.beta =pbeta; 
results.rho = prho; 
bhat = results.beta;
results.sig2_v = psig2_v; 
results.Gamma=pGamma;
results.sig2_ep=psig2_ep;
results.coeffcorr=pcoeffcorr;

%%%%%%%%%%%%%%%%%%%%%%%%%%

yhat = (speye(n) - prho*W)\(x1*bhat);
results.yhat = yhat;
results.resid = y - yhat;


parm= [results.rho
       results.beta
       results.Gamma
       results.sig2_ep
       results.sig2_v
       results.coeffcorr
       ];      



results.lik = f2_sar_Wendo(parm,detval,n,Zn,x2,y,W,x1);


% asymptotic t-stats using numerical hessian
t0 = clock;
% just computes the diagonal
f2_sar_Wendo_oc=@(x) f2_sar_Wendo(x,detval,n,Zn,x2,y,W,x1);
dhessn = hessian(f2_sar_Wendo_oc,parm);
%dhessn = hessnum;
hessi_ini = invpd(-dhessn);

%tousd les coeff estim�: varcov et estimateurs
results.Allvarcov=hessi_ini;
results.Allparam=parm;




%On ne prend que les estimateur dont on a besoin: rho et beta
hessi1=hessi_ini(1:size(x1,2)+1,1:size(x1,2)+1);
%On r�arrange en mettant rho � la fin pour avoir les mm notations que le code de Lesage
hessi2=[hessi1(:,2:end) hessi1(:,1)];
hessi3=[hessi2(2:end,:); hessi2(1,:)];
hessi=hessi3;



if any(diag(hessi) < 0)
    fprintf(1,'sar: negative variances from numerical hessian \n');
    fprintf(1,'sar: t-statistics may be inaccurate \n');
end;




tvar = abs(diag(hessi));




tmp = [results.beta
    results.rho
       ];
   


results.tstat = tmp./sqrt(tvar);
time3 = etime(clock,t0);
results.time3 = time3;

% end of t-stat calculations


%% Impacts

% do effects estimates

t0 = clock; 

% pre-calculate traces for the x-impacts calculations
uiter=50;
maxorderu=100;
nobs = n;
rv=randn(nobs,uiter);
tracew=zeros(maxorderu,1);
wjjju=rv;
for jjj=1:maxorderu
    wjjju=W*wjjju;
    tracew(jjj)=mean(mean(rv.*wjjju));
    
end

traces=[tracew];
traces(1)=0;
traces(2)=sum(sum(W'.*W))/nobs;
trs=[1;traces];
ntrs=length(trs);
trbig=trs';

% cheat here to fix the numerical hessian if it sucks
% Use MCMC to get good results
sigmat = hessi - diag(diag(hessi)) + diag(diag(abs(hessi)));
sigmatt = sigmat; %On en prend que rho et beta pour simulation monte carlo
[R,posdef] = chol(sigmatt);



tmp = [results.beta
       results.rho];
   
bdraws = matadd(norm_rndmat(sigmatt,ndraw),tmp);
draws = bdraws';

psave = draws(:,end);
ind = find(psave > 1); % find bad rho draws
psave(ind,1) = 0.99; % replace them with 0.99


bsave = draws(:,1:end-1);

        if cflag == 1
        bdraws = bsave(:,2:end);
        elseif cflag == 0
        bdraws = bsave;
        end; 
        pdraws = psave;
        
        ree = 0:1:ntrs-1;

        rmat = zeros(1,ntrs);
        total = zeros(ndraw,p,ntrs);
        direct = zeros(ndraw,p,ntrs);
        indirect = zeros(ndraw,p,ntrs);

        
for i=1:ndraw;
    rmat = pdraws(i,1).^ree;
    for j=1:p;
            beta = [bdraws(i,j)];
            total(i,j,:) = beta(1,1)*rmat;
    direct(i,j,:) = (beta*trbig).*rmat;
    indirect(i,j,:) = total(i,j,:) - direct(i,j,:);
    end;

end;

time5 = etime(clock,t0);
results.time5 = time5;




ym = y - mean(y);       % r-squared, rbar-squared
rsqr1 = results.resid'*results.resid;
rsqr2 = ym'*ym;
results.rsqr = 1.0-rsqr1/rsqr2;   % r-squared
rsqr1 = rsqr1/(n-nvar);
rsqr2 = rsqr2/(n-1.0);

% return stuff
results.meth = 'sar_Wendo';
results.y = y;      
results.total = total;
results.direct = direct;
results.indirect = indirect;
results.nobs = n; 
results.nvar = nvar;
results.rmax = rmax;      
results.rmin = rmin;
results.lflag = ldetflag;
results.order = order;
results.miter = miter;
results.rbar = 1 - (rsqr1/rsqr2); % rbar-squared
results.time = etime(clock,timet);
results.lndet = detval;






function llike = f2_sar_Wendo(theta,detval,n,Zn,x2,y,W,x1)
% PURPOSE: evaluates log-likelihood -- given ML estimates
%  spatial autoregressive model using sparse matrix algorithms
% ---------------------------------------------------
%  USAGE:llike = f2_sar(parm,y,X,W,ldet)
%  where: parm = vector of maximum likelihood parameters
%                parm(1:k-2,1) = b, parm(k-1,1) = rho, parm(k,1) = sige
%         y    = dependent variable vector (n x 1)
%         X    = explanatory variables matrix (n x k)
%         W    = spatial weight matrix
%         ldet = matrix with [rho log determinant] values
%                computed in sar.m using one of Kelley Pace's routines  
% ---------------------------------------------------
%  RETURNS: a  scalar equal to minus the log-likelihood
%           function value at the ML parameters
%  --------------------------------------------------
%  NOTE: this is really two functions depending
%        on nargin = 4 or nargin = 5 (see the function)
% ---------------------------------------------------
%  SEE ALSO: sar, f2_far, f2_sac, f2_sem
% ---------------------------------------------------

% written by: James P. LeSage 1/2000
% University of Toledo
% Department of Economics
% Toledo, OH 43606
% jlesage@spatial.econometrics.com

rho=theta(1,1);
beta=theta(2:size(x1,2)+1,1);
Gamma=theta(size(x1,2)+2:size(x1,2)+1+size(x2,2),1); %hypoth�se que p2=1
sig2_ep=theta(size(x1,2)+2+size(x2,2),1); %hypoth�se que p2=1
sig2_v=theta(size(x1,2)+2+size(x2,2)+1,1);
coeffcorr=theta(size(x1,2)+2+size(x2,2)+2,1); %hypoth�se que p2=1



gsize = detval(2,1) - detval(1,1);
i1 = find(detval(:,1) <= rho + gsize);
i2 = find(detval(:,1) <= rho - gsize);
i1 = max(i1);
i2 = max(i2);
index = round((i1+i2)/2);
if isempty(index)
index = 1;
end;
detm = detval(index,2);


ep=Zn-x2*Gamma;
eta=(eye(n) - rho*W)*y-x1*beta-ep*coeffcorr*sqrt(sig2_v)/sqrt(sig2_ep);


llike = -(n)*log(2*pi) -(n/2)*log(sig2_ep*sig2_v*(1-coeffcorr^2)) + detm - 1/(2*sig2_ep)*(ep')*ep - 1/(2*sig2_v*(1-coeffcorr^2))*(eta')*eta;


function [rmin,rmax,convg,maxit,detval,ldetflag,eflag,order,iter,options,ndraw] = sar_parse(info)
% PURPOSE: parses input arguments for sar model
% ---------------------------------------------------
%  USAGE: [rmin,rmax,convg,maxit,detval,ldetflag,eflag,order,iter,options,ndraw] = sar_parse(info)
% where info contains the structure variable with inputs 
% and the outputs are either user-inputs or default values
% ---------------------------------------------------

% set defaults
options = zeros(1,18); % optimization options for fminbnd
options(1) = 0; 
options(2) = 1.e-6; 
options(14) = 500;

eflag = 0;     % default to not computing eigenvalues
ldetflag = 1;  % default to 1999 Pace and Barry MC determinant approx
order = 50;    % there are parameters used by the MC det approx
iter = 30;     % defaults based on Pace and Barry recommendation
rmin = -1;     % use -1,1 rho interval as default
rmax = 1;
detval = 0;    % just a flag
convg = 0.0001;
maxit = 500;
ndraw = 1000;


fields = fieldnames(info);
nf = length(fields);
if nf > 0
    
 for i=1:nf
    if strcmp(fields{i},'rmin')
        rmin = info.rmin;  eflag = 0;
    elseif strcmp(fields{i},'rmax')
        rmax = info.rmax; eflag = 0;
    elseif strcmp(fields{i},'convg')
        options(2) = info.convg;
    elseif strcmp(fields{i},'maxit')
        options(14) = info.maxit;  
    elseif strcmp(fields{i},'lndet')
    detval = info.lndet;
    ldetflag = -1;
    eflag = 0;
    rmin = detval(1,1);
    nr = length(detval);
    rmax = detval(nr,1);
    elseif strcmp(fields{i},'lflag')
        tst = info.lflag;
        if tst == 0,
        ldetflag = 0; % compute full lndet, no approximation
        elseif tst == 1,
        ldetflag = 1; % use Pace-Barry approximation
        elseif tst == 2,
        ldetflag = 2; % use spline interpolation approximation
        else
        error('sar: unrecognizable lflag value on input');
        end;
    elseif strcmp(fields{i},'order')
        order = info.order;  
    elseif strcmp(fields{i},'eig')
        eflag = info.eig;  
    elseif strcmp(fields{i},'iter')
        iter = info.iter; 
     elseif strcmp(fields{i},'ndraw')
        ndraw = info.ndraw; 
    end;
 end;
 
else % the user has input a blank info structure
      % so we use the defaults
end; 

function [rmin,rmax,time2] = sar_eigs(eflag,W,rmin,rmax,n)
% PURPOSE: compute the eigenvalues for the weight matrix
% ---------------------------------------------------
%  USAGE: [rmin,rmax,time2] = far_eigs(eflag,W,rmin,rmax,W)
% where eflag is an input flag, W is the weight matrix
%       rmin,rmax may be used as default outputs
% and the outputs are either user-inputs or default values
% ---------------------------------------------------


if eflag == 1 % do eigenvalue calculations
t0 = clock;
opt.tol = 1e-3; opt.disp = 0;
lambda = eigs(sparse(W),speye(n),1,'SR',opt);  
rmin = real(1/lambda);   
rmax = 1.0;
time2 = etime(clock,t0);
else % use rmin,rmax arguments from input or defaults -1,1
time2 = 0;
end;


function [detval,time1] = sar_lndet(ldetflag,W,rmin,rmax,detval,order,iter)
% PURPOSE: compute the log determinant |I_n - rho*W|
% using the user-selected (or default) method
% ---------------------------------------------------
%  USAGE: detval = far_lndet(lflag,W,rmin,rmax)
% where eflag,rmin,rmax,W contains input flags 
% and the outputs are either user-inputs or default values
% ---------------------------------------------------


% do lndet approximation calculations if needed
if ldetflag == 0 % no approximation
t0 = clock;    
out = lndetfull(W,rmin,rmax);
time1 = etime(clock,t0);
tt=rmin:.001:rmax; % interpolate a finer grid
outi = interp1(out.rho,out.lndet,tt','spline');
detval = [tt' outi];
    
elseif ldetflag == 1 % use Pace and Barry, 1999 MC approximation

t0 = clock;    
out = lndetmc(order,iter,W,rmin,rmax);
time1 = etime(clock,t0);
results.limit = [out.rho out.lo95 out.lndet out.up95];
tt=rmin:.001:rmax; % interpolate a finer grid
outi = interp1(out.rho,out.lndet,tt','spline');
detval = [tt' outi];

elseif ldetflag == 2 % use Pace and Barry, 1998 spline interpolation

t0 = clock;
out = lndetint(W,rmin,rmax);
time1 = etime(clock,t0);
tt=rmin:.001:rmax; % interpolate a finer grid
outi = interp1(out.rho,out.lndet,tt','spline');
detval = [tt' outi];

elseif ldetflag == -1 % the user fed down a detval matrix
    time1 = 0;
        % check to see if this is right
        if detval == 0
            error('sar: wrong lndet input argument');
        end;
        [n1,n2] = size(detval);
        if n2 ~= 2
            error('sar: wrong sized lndet input argument');
        elseif n1 == 1
            error('sar: wrong sized lndet input argument');
        end;          
end;



function H = hessian(f,x,varargin)
% PURPOSE: Computes finite difference Hessian
% -------------------------------------------------------
% Usage:  H = hessian(func,x,varargin)
% Where: func = function name, fval = func(x,varargin)
%           x = vector of parameters (n x 1)
%    varargin = optional arguments passed to the function
% -------------------------------------------------------
% RETURNS:
%           H = finite differnce hessian
% -------------------------------------------------------

% Code from:
% COMPECON toolbox [www4.ncsu.edu/~pfackler]
% documentation modified to fit the format of the Ecoometrics Toolbox
% by James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jlesage@spatial-econometrics.com

eps = 1e-8;

n = size(x,1);
fx = feval(f,x,varargin{:});
 
% Compute the stepsize (h)
h = eps.^(1/3)*max(abs(x),1e-2);
xh = x+h;
h = xh-x;    
ee = sparse(1:n,1:n,h,n,n);
 
% Compute forward step 
g = zeros(n,1);
for i=1:n
  g(i) = feval(f,x+ee(:,i),varargin{:});
end
   
H=h*h';
% Compute "double" forward step 
for i=1:n
for j=i:n
  H(i,j) = (feval(f,x+ee(:,i)+ee(:,j),varargin{:})-g(i)-g(j)+fx)/H(i,j);
  H(j,i) = H(i,j);
end
end

function y = norm_rndmat(sig,ndraw)
% PURPOSE: random multivariate random vector based on
%          var-cov matrix sig
%---------------------------------------------------
% USAGE:   y = norm_rnd(sig)
% where:   sig = a square-symmetric covariance matrix 
% NOTE: for mean b, var-cov sig use: b +  norm_rnd(sig) 
%---------------------------------------------------      
% RETURNS: y = random vector normal draw mean 0, var-cov(sig)
%---------------------------------------------------

% by 
% James P. LeSage, last updated 3/2010
% Dept of Finance & Economics
% Texas State University-San Marcos
% 601 University Drive
% San Marcos, TX 78666
% jlesage@spatial-econometrics.com

if nargin ~= 2
error('Wrong # of arguments to norm_rnd');
end;

h = chol(sig);
[nrow, ncol] = size(sig);
rv = randn(nrow,ndraw);

y = h'*rv;



