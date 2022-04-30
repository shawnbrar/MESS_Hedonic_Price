function results = sar_gv(y,x,W,ndraw,nomit,prior)
% PURPOSE: Bayesian estimates of the spatial autoregressive model
% THIS FUNCTION: estimates the heteroscedasticity parameter r
%                and returns results in results.rdraw for posterior
%                inference
%          y = rho*W*y + XB + e, e = N(0,sige*V), V = diag(v1,v2,...vn) 
%          r/vi = ID chi(r)/r
%          r = Gamma(delta,2)
%          B = N(c,T), 
%          1/sige = Gamma(nu,d0), 
%          rho = Uniform(rmin,rmax), or rho = beta(a1,a2); 
%-------------------------------------------------------------
% USAGE: results = sar_gv(y,x,W,ndraw,nomit,prior)
% where: y = dependent variable vector (nobs x 1)
%        x = independent variables matrix (nobs x nvar)
%        W = spatial weight matrix (standardized, row-sums = 1)
%    ndraw = # of draws
%    nomit = # of initial draws omitted for burn-in            
%    prior = a structure variable with:
%            prior.beta  = prior means for beta,   c above (default 0)
%            priov.bcov  = prior beta covariance , T above (default 1e+12)
%            prior.nu    = informative Gamma(nu,d0) prior on sige
%            prior.d0    = default: nu=0,d0=0 (diffuse prior)
%            prior.delta = default: delta = 20
%            prior.a1    = parameter for beta(a1,a2) prior on rho see: 'help beta_prior'
%            prior.a2    = (default = 1.0, a uniform prior on rmin,rmax) 
%            prior.eig   = 0 for default rmin = -1,rmax = +1, 1 for eigenvalue calculation of these
%            prior.rmin  = (optional) min rho used in sampling (default = -1)
%            prior.rmax  = (optional) max rho used in sampling (default = 1)  
%            prior.lflag = 0 for full lndet computation (default = 1, fastest)
%                        = 1 for MC approx (fast for large problems)
%                        = 2 for Spline approx (medium speed)
%            prior.order = order to use with prior.lflag = 1 option (default = 50)
%            prior.iter  = iters to use with prior.lflag = 1 option (default = 30) 
%            prior.lndet = a matrix returned by sar, sar_g, sarp_g, etc.
%                          containing log-determinant information to save time
%            prior.logm  = 0 for no log marginal calculation, = 1 for log marginal (default = 1)
%            prior.sflag = 1 if this function is being called by sdm_g(), default to not used
%-------------------------------------------------------------
% RETURNS:  a structure:
%          results.meth     = 'sar_g'
%          results.beta     = posterior mean of bhat based on draws
%          results.rho      = posterior mean of rho based on draws
%          results.sige     = posterior mean of sige based on draws
%          results.sigma    = posterior mean of sige based on (e'*e)/(n-k)
%          results.bdraw    = bhat draws (ndraw-nomit x nvar)
%          results.pdraw    = rho  draws (ndraw-nomit x 1)
%          results.rdraw    = r draws (ndraw-nomit x 1) 
%          results.sdraw    = sige draws (ndraw-nomit x 1)
%          results.total    = a 3-d matrix (ndraw,nvars-1,ntrs) total x-impacts
%          results.direct   = a 3-d matrix (ndraw,nvars-1,ntrs) direct x-impacts
%          results.indirect = a 3-d matrix (ndraw,nvars-1,ntrs) indirect x-impacts
%                             ntrs defaults to 101 trace terms
%          results.vmean    = mean of vi draws (nobs x 1) 
%          results.delta  = value of delta (from input)
%          results.bmean  = b prior means, prior.beta from input
%          results.bstd   = b prior std deviations sqrt(diag(prior.bcov))
%          results.nobs   = # of observations
%          results.nvar   = # of variables in x-matrix
%          results.ndraw  = # of draws
%          results.nomit  = # of initial draws omitted
%          results.y      = y-vector from input (nobs x 1)
%          results.yhat   = mean of posterior predicted (nobs x 1)
%          results.resid  = residuals, based on posterior means
%          results.rsqr   = r-squared based on posterior means
%          results.rbar   = adjusted r-squared
%          results.nu     = nu prior parameter
%          results.d0     = d0 prior parameter
%          results.a1     = a1 parameter for beta prior on rho from input, or default value
%          results.a2     = a2 parameter for beta prior on rho from input, or default value
%          results.time1  = time for eigenvalue calculation
%          results.time2  = time for log determinant calcluation
%          results.time3  = time for sampling
%          results.time   = total time taken  
%          results.rmax   = 1/max eigenvalue of W (or rmax if input)
%          results.rmin   = 1/min eigenvalue of W (or rmin if input)          
%          results.tflag  = 'plevel' (default) for printing p-levels
%                         = 'tstat' for printing bogus t-statistics 
%          results.lflag  = lflag from input
%          results.iter   = prior.iter option from input
%          results.order  = prior.order option from input
%          results.limit  = matrix of [rho lower95,logdet approx, upper95] 
%                           intervals for the case of lflag = 1
%          results.lndet = a matrix containing log-determinant information
%                          (for use in later function calls to save time)
%          results.mlike = log marginal likelihood (a vector ranging over
%                          rho values that can be integrated for model comparison)
%          results.cflag  = 1 for intercept term, 0 for no intercept term
% --------------------------------------------------------------
% NOTES: - use either improper prior.rval 
%          or informative Gamma prior.m, prior.k, not both of them
% - for n < 1000 you should use lflag = 0 to get exact results  
% - use a1 = 1.0 and a2 = 1.0 for uniform prior on rho
% - results.mlike can be used for model comparison (see model_compare.m for a demonstration)
% --------------------------------------------------------------
% SEE ALSO: (sar_gd, sar_gd2 demos) prt
% --------------------------------------------------------------
% REFERENCES: LeSage and Pace (2009) Chapter 5 on Bayesian estimation 
%             of spatial regression models.
% For lndet information see: Chapter 4 
%----------------------------------------------------------------

% written by:
% James P. LeSage, last updated 3/2009
% Dept of Finance & Economics
% Texas State University-San Marcos
% 601 University Drive
% San Marcos, TX 78666
% jlesage@spatial-econometrics.com


timet = clock;

% error checking on inputs
[n junk] = size(y);
[n1 k] = size(x);
[n3 n4] = size(W);
time1 = 0;
time2 = 0;
time3 = 0;

nobsa = n;

results.nobs  = n;
results.nvar  = k;
results.y = y; 

n = length(y);
if sum(x(:,1)) ~= n
tst = sum(x); % we may have no intercept term
ind = find(tst == n); % we do have an intercept term
if length(ind) > 0
error('sar_g: intercept term must be in first column of the x-matrix');
elseif length(ind) == 0 % case of no intercept term
cflag = 0;
p = size(x,2);
end;
elseif sum(x(:,1)) == n % we have an intercept in the right place
cflag = 1;
p = size(x,2)-1;
end;

if nargin == 5
    prior.lflag = 1;
end;
  
[nu,d0,rval,mm,kk,rho,sige,rmin,rmax,detval,ldetflag,sflag, ...
eflag,order,iter,c,T,inform_flag,a1,a2,logmflag,delta] = sar_parse(prior,k);

results.delta = delta;

if sflag == 1
p = (k-1)/2;
if n1 ~= n
error('sdm_g: x-matrix contains wrong # of observations');
elseif n3 ~= n4
error('sdm_g: W matrix is not square');
elseif n3~= n
error('sdm_g: W matrix is not the same size at y,x');
end;

% error checking on prior information inputs
[checkk,junk] = size(c);
if checkk ~= k
error('sdm_g: prior means are wrong');
elseif junk ~= 1
error('sdm_g: prior means are wrong');
end;

[checkk junk] = size(T);
if checkk ~= k
error('sdm_g: prior bcov is wrong');
elseif junk ~= k
error('sdm_g: prior bcov is wrong');
end;


elseif sflag == 0

p = k-1;

if n1 ~= n
error('sar_g: x-matrix contains wrong # of observations');
elseif n3 ~= n4
error('sar_g: W matrix is not square');
elseif n3~= n
error('sar_g: W matrix is not the same size at y,x');
end;

% error checking on prior information inputs
[checkk,junk] = size(c);
if checkk ~= k
error('sar_g: prior means are wrong');
elseif junk ~= 1
error('sar_g: prior means are wrong');
end;

[checkk junk] = size(T);
if checkk ~= k
error('sar_g: prior bcov is wrong');
elseif junk ~= k
error('sar_g: prior bcov is wrong');
end;

end



results.order = order;
results.iter = iter;

timet = clock; % start the timer

[rmin,rmax,time1] = sar_eigs(eflag,W,rmin,rmax,n);

[detval,time2] = sar_lndet(ldetflag,W,rmin,rmax,detval,order,iter);

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


trs=[1;traces];
ntrs=length(trs);
strs=(0:ntrs)';
op=ones(p,1);
trbig=trs(:,op)';
maxorder = 1;
ereject = 0;

% storage for draws
          bsave = zeros(ndraw-nomit,k);
          rsave = zeros(ndraw-nomit,1);
          psave = zeros(ndraw-nomit,1);
          ssave = zeros(ndraw-nomit,1);
          vmean = zeros(n,1);
          acc_rate = zeros(ndraw,1);
          ccsave = zeros(ndraw,1);

          total   =zeros(ndraw-nomit,p,ntrs);
          direct  =zeros(ndraw-nomit,p,ntrs);
          indirect=zeros(ndraw-nomit,p,ntrs);

% ====== initializations
% compute this stuff once to save time
TI = inv(T);
TIc = TI*c;
iter = 1;

in = ones(n,1);
V = in;
vi = in;
Wy = sparse(W)*y;

%Prior for degrees of freedom is exponential with mean vl0
vl0=delta;
cc = 5;
pswitch = 0;
acc = 0;


    
if sflag == 0
hwait = waitbar(0,'sar: MCMC sampling ...');
elseif sflag == 1
hwait = waitbar(0,'sdm: MCMC sampling ...');
end;

t0 = clock;                  
iter = 1;
          while (iter <= ndraw); % start sampling;
                  
          % update beta   
          xs = matmul(x,sqrt(V));
          ys = sqrt(V).*y;
          Wys = sqrt(V).*Wy;
          AI = inv(xs'*xs + sige*TI);         
          yss = ys - rho*Wys;          
          xpy = xs'*yss;
          b = xs'*yss + sige*TIc;
          b0 = AI*b;
          bhat = norm_rnd(sige*AI) + b0;  
          xb = xs*bhat;
                    
          % update sige
          nu1 = nobsa + 2*nu; 
          e = (yss - xb);
          d1 = 2*d0 + e'*e;
          chi = chis_rnd(1,nu1);
          sige = d1/chi;
          
          % update vi
          ev = ys - rho*Wys - xs*bhat; 
          dof = vdraw + 1;
          error2 = ev.*ev;
          tmp = (1/sige)*error2 + vdraw;
          chiv = chis_rnd(n,dof);   
          V = chiv./tmp;

          
      % we use griddy Gibbs to perform rho-draw
          b0 = (xs'*xs + sige*TI )\(xs'*ys + sige*TIc);
          bd = (xs'*xs + sige*TI)\(xs'*Wys + sige*TIc);
          e0 = ys - xs*b0;
          ed = Wys - xs*bd;
          epe0 = e0'*e0;
          eped = ed'*ed;
          epe0d = ed'*e0;
          logdetx = log(det(xs'*xs + TI));
          rho = draw_rho(detval,epe0,eped,epe0d,n,k,rho,a1,a2,logdetx);
          
               %Random walk Metropolis step for dof
    temp = -log(V) + V;
    nu = 1/vl0 + .5*sum(temp);
    
     vlcan= vdraw +  cc*randn(1,1);
     if vlcan>0
        lpostcan = .5*n*vlcan*log(.5*vlcan) -n*gammaln(.5*vlcan)...
        -nu*vlcan;
        lpostdraw = .5*n*vdraw*log(.5*vdraw) -n*gammaln(.5*vdraw)...
        -nu*vdraw;
        accprob = exp(lpostcan-lpostdraw);
     else
        accprob=0;
     end
     

%accept candidate draw with log prob = laccprob, else keep old draw
   if  rand<accprob
       vdraw=vlcan;
       pswitch=pswitch+1;
       acc = acc + 1;
   end    

      acc_rate(iter,1) = acc/iter;
      % update cc based on std of rho draws
       if acc_rate(iter,1) < 0.4
       cc = cc/1.1;
       ccsave(iter,1) = cc;
       end;
       if acc_rate(iter,1) > 0.6
       cc = cc*1.1;
       ccsave(iter,1) = cc;
       end;
   
   

        % calculate effects estimates
        if sflag == 1
        maxorder1=maxorder+1;
        [totale,directe,indirecte] = calc_sdm_effects(maxorder1,trs,strs,ntrs,rho,bhat,p);
        elseif sflag == 0
        [totale,directe,indirecte] = calc_sar_effects(maxorder,trs,strs,ntrs,rho,bhat,p);
        end;


          
    if iter > nomit % if we are past burn-in, save the draws
        bsave(iter-nomit,1:k) = bhat';
        ssave(iter-nomit,1) = sige;
        psave(iter-nomit,1) = rho;
        rsave(iter-nomit,1) = vdraw;

        vmean = vmean + in./V;
        total(iter-nomit,:,:)=totale; % a p by ntraces matrix
        direct(iter-nomit,:,:)=directe; % a p by ntraces matrix
        indirect(iter-nomit,:,:)=indirecte; % a p by ntraces matrix         
    

    if mm~= 0
        rsave(iter-nomit,1) = rval;
    end;         
    end;
                    
iter = iter + 1; 
waitbar(iter/ndraw);         
end; % end of sampling loop
close(hwait);

time3 = etime(clock,t0);

% compute posterior means and log marginal likelihood for return arguments
bmean = mean(bsave);
beta = bmean';
rho = mean(psave);
sige = mean(ssave);
results.sige = sige;
vmean = vmean/(ndraw-nomit);
V = in./vmean;
[nobs,nvar] = size(x);
if logmflag == 1
          xs = matmul(x,sqrt(V));
          ys = sqrt(V).*y;
          Wys = W*ys;
          AI = inv(xs'*xs + sige*TI);
          b0 = AI*(xs'*ys + sige*TIc);
          bd = AI*(xs'*Wys + sige*TIc);
          e0 = ys - xs*b0;
          ed = Wys - xs*bd;
          epe0 = e0'*e0;
          eped = ed'*ed;
          epe0d = ed'*e0;
 logdetx = log(det(xs'*xs + sige*TI));
 mlike = sar_marginal(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,logdetx,a1,a2);
 yhat = (speye(nobs) - rho*W)\(x*beta);
 e = y - yhat;
elseif logmflag == 0
 xs = matmul(x,sqrt(V));
 ys = sqrt(V).*y;
 yhat = (speye(nobs) - rho*W)\(x*beta);
 e = y - yhat;
end;
% compute R-squared
epe = e'*e;
sige = epe/(n-k);
results.sigma = sige;
ym = y - mean(y);
rsqr1 = epe;
rsqr2 = ym'*ym;
results.rsqr = 1- rsqr1/rsqr2; % r-squared
rsqr1 = rsqr1/(nobs-nvar);
rsqr2 = rsqr2/(nobs-1.0);
results.rbar = 1 - (rsqr1/rsqr2); % rbar-squared

time = etime(clock,timet);


results.meth  = 'sar_g';
results.total = total;
results.direct = direct;
results.indirect = indirect;
results.beta = beta;
results.rho = rho;
results.bdraw = bsave;
results.pdraw = psave;
results.sdraw = ssave;
if logmflag == 1
results.mlike = mlike;
end;
results.vmean = vmean;
results.yhat  = yhat;
results.resid = e;
results.bmean = c;
results.bstd  = sqrt(diag(T));
results.ndraw = ndraw;
results.nomit = nomit;
results.time  = time;
results.time1 = time1;
results.time2 = time2;
results.time3 = time3;
results.nu = nu;
results.d0 = d0;
results.a1 = a1;
results.a2 = a2;
results.tflag = 'plevel';
results.rmax = rmax; 
results.rmin = rmin;
results.lflag = ldetflag;
results.lndet = detval;
results.novi  = novi_flag;
results.priorb = inform_flag;

results.rdraw = rsave;


% =========================================================================
% support functions below
% =========================================================================

function rho = draw_rho(detval,epe0,eped,epe0d,n,k,rho,a1,a2,logdetx)

nmk = (n-k)/2;
nrho = length(detval(:,1));
iota = ones(nrho,1);

z = epe0*iota - 2*detval(:,1)*epe0d + detval(:,1).*detval(:,1)*eped;
z = -nmk*log(z);
%C = gammaln(nmk)*iota -nmk*log(2*pi)*iota - 0.5*logdetx*iota;
den =  detval(:,2) + z;

bprior = beta_prior(detval(:,1),a1,a2);
den = den + log(bprior);
n = length(den);
y = detval(:,1);
adj = max(den);
den = den - adj;
x = exp(den);
% trapezoid rule
isum = sum((y(2:n,1) + y(1:n-1,1)).*(x(2:n,1) - x(1:n-1,1))/2);
z = abs(x/isum);
den = cumsum(z);

rnd = unif_rnd(1,0,1)*sum(z);
ind = find(den <= rnd);
idraw = max(ind);
if (idraw > 0 & idraw < nrho)
rho = detval(idraw,1);
end;

% To see how this works, uncomment the following lines
% plot(detval(:,1),den/1000,'-');
% line([detval(idraw,1) detval(idraw,1)],[0 den(idraw,1)/1000]);
% hold on;
% line([detval(idraw,1) 0],[den(idraw,1)/1000 den(idraw,1)/1000]);
% drawnow;
% pause;


function [nu,d0,rval,mm,kk,rho,sige,rmin,rmax,detval,ldetflag,sflag,eflag,order,iter,c,T,inform_flag,a1,a2,logmflag,delta] = sar_parse(prior,k)
% PURPOSE: parses input arguments for sar_g models
% ---------------------------------------------------
%  USAGE: [nu,d0,rval,mm,kk,rho,sige,rmin,rmax,detval, ...
%         ldetflag,eflag,mflag,order,iter,novi_flag,c,T,inform_flag,a1,a2,logmflag = 
%                           sar_parse(prior,k)
% where info contains the structure variable with inputs 
% and the outputs are either user-inputs or default values
% ---------------------------------------------------

% set defaults

logmflag = 1;  % default to compute log-marginal
eflag = 0;     % default to not computing eigenvalues
ldetflag = 1;  % default to 1999 Pace and Barry MC determinant approx
mflag = 1;     % default to compute log marginal likelihood
order = 50;    % there are parameters used by the MC det approx
iter = 30;     % defaults based on Pace and Barry recommendation
rmin = -1;     % use -1,1 rho interval as default
rmax = 1;
detval = 0;    % just a flag
rho = 0.5;
sige = 1.0;
rval = 4;
mm = 0;
kk = 0;
nu = 0;
d0 = 0;
a1 = 1.0;
a2 = 1.0;
c = zeros(k,1);   % diffuse prior for beta
T = eye(k)*1e+12;
prior_beta = 0;   % flag for diffuse prior on beta
inform_flag = 0;
sflag = 0; % default to sar model
delta = 20;

fields = fieldnames(prior);
nf = length(fields);
if nf > 0
 for i=1:nf
    if strcmp(fields{i},'nu')
        nu = prior.nu;
    elseif strcmp(fields{i},'d0')
        d0 = prior.d0;  
    elseif strcmp(fields{i},'delta')
        delta = prior.delta; 
    elseif strcmp(fields{i},'logm')
       logmflag = prior.logm; 
    elseif strcmp(fields{i},'a1')
       a1 = prior.a1; 
    elseif strcmp(fields{i},'a2')
       a2 = prior.a2; 
    elseif strcmp(fields{i},'beta')
        c = prior.beta; inform_flag = 1; % flag for informative prior on beta
    elseif strcmp(fields{i},'bcov')
        T = prior.bcov; inform_flag = 1;
    elseif strcmp(fields{i},'rmin')
        rmin = prior.rmin; eflag = 0;
    elseif strcmp(fields{i},'rmax')
        rmax = prior.rmax;  eflag = 0;
    elseif strcmp(fields{i},'lndet')
    detval = prior.lndet;
    ldetflag = -1;
    eflag = 0;
    rmin = detval(1,1);
    nr = length(detval);
    rmax = detval(nr,1);
    elseif strcmp(fields{i},'lflag')
        tst = prior.lflag;
        if tst == 0,
        ldetflag = 0; 
        elseif tst == 1,
        ldetflag = 1; 
        elseif tst == 2,
        ldetflag = 2; 
        else
        error('sar_g: unrecognizable lflag value on input');
        end;
    elseif strcmp(fields{i},'order')
        order = prior.order;  
    elseif strcmp(fields{i},'iter')
        iter = prior.iter; 
    elseif strcmp(fields{i},'dflag')
        metflag = prior.dflag;
    elseif strcmp(fields{i},'eig')
        eflag = prior.eig;
    elseif strcmp(fields{i},'sflag')
        sflag = prior.sflag;
    end;
 end;

 
else, % the user has input a blank info structure
      % so we use the defaults
end; 

function [rmin,rmax,time2] = sar_eigs(eflag,W,rmin,rmax,n);
% PURPOSE: compute the eigenvalues for the weight matrix
% ---------------------------------------------------
%  USAGE: [rmin,rmax,time2] = far_eigs(eflag,W,rmin,rmax,W)
% where eflag is an input flag, W is the weight matrix
%       rmin,rmax may be used as default outputs
% and the outputs are either user-inputs or default values
% ---------------------------------------------------


if eflag == 1 % compute eigenvalues
t0 = clock;
opt.tol = 1e-3; opt.disp = 0;
lambda = eigs(sparse(W),speye(n),1,'SR',opt);  
rmin = 1/real(lambda);   
rmax = 1;
time2 = etime(clock,t0);
else
time2 = 0;
end;


function [detval,time1] = sar_lndet(ldetflag,W,rmin,rmax,detval,order,iter);
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
            error('sar_g: wrong lndet input argument');
        end;
        [n1,n2] = size(detval);
        if n2 ~= 2
            error('sar_g: wrong sized lndet input argument');
        elseif n1 == 1
            error('sar_g: wrong sized lndet input argument');
        end;          
end;

function  out = sar_marginal(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,logdetx,a1,a2)
% PURPOSE: returns a vector of the log-marginal over a grid of rho-values
% -------------------------------------------------------------------------
% USAGE: out = sar_marginal(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,logdetx,a1,a2)
% where:       detval = an ngrid x 2 matrix with rho-values and lndet values
%                  e0 = y - x*b0;
%                 ed = Wy - x*bd;
%               epe0 = e0'*e0;
%               eped = ed'*ed;
%              epe0d = ed'*e0;
%               nobs = # of observations
%               nvar = # of explanatory variables
%            logdetx = log(det(x'*x))
%                 a1 = parameter for beta prior on rho
%                 a2 = parameter for beta prior on rho
% -------------------------------------------------------------------------
% RETURNS: out = a structure variable
%        out = log marginal, a vector the length of detval
% -------------------------------------------------------------------------

% written by:
% James P. LeSage, 7/2003
% Dept of Finance & Economics
% Texas State University-San Marcos
% 601 University Drive
% San Marcos, TX 78666
% jlesage@spatial-econometrics.com

n = length(detval);
nmk = (nobs-nvar)/2;
% C is a constant of integration that can vary with nvars, so for model
% comparisions involving different nvars we need to include this
bprior = beta_prior(detval(:,1),a1,a2);
C = log(bprior) + gammaln(nmk) - nmk*log(2*pi) - 0.5*logdetx;
iota = ones(n,1);
z = epe0*iota - 2*detval(:,1)*epe0d + detval(:,1).*detval(:,1)*eped;
den = C + detval(:,2) - nmk*log(z);
out = real(den);

function  out = sar_marginal2(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,a1,a2,c,TI,xs,ys,sige,W)
% PURPOSE: returns a vector of the log-marginal over a grid of rho-values
%          for the case of an informative prior on beta
% -------------------------------------------------------------------------
% USAGE: out = sar_marginal(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,logdetx,a1,a2)
% where:       detval = an ngrid x 2 matrix with rho-values and lndet values
%                  e0 = y - x*b0;
%                 ed = Wy - x*bd;
%               epe0 = e0'*e0;
%               eped = ed'*ed;
%              epe0d = ed'*e0;
%               nobs = # of observations
%               nvar = # of explanatory variables
%            logdetx = log(det(x'*x))
%                 a1 = parameter for beta prior on rho
%                 a2 = parameter for beta prior on rho
%                 c = prior mean for beta
%                TI = prior var-cov for beta
%                xs = x*sqrt(V) or x if homoscedastic model
%                ys = y*sqrt(V) or y is homoscedastic model
% -------------------------------------------------------------------------
% RETURNS: out = a structure variable
%        out = log marginal, a vector the length of detval
% -------------------------------------------------------------------------

% written by:
% James P. LeSage, 7/2003
% Dept of Finance & Economics
% Texas State University-San Marcos
% 601 University Drive
% San Marcos, TX 78666
% jlesage@spatial-econometrics.com

n = length(detval);
nmk = (nobs-nvar)/2;
% C is a constant of integration that can vary with nvars, so for model
% comparisions involving different nvars we need to include this
bprior = beta_prior(detval(:,1),a1,a2);
C = log(bprior) + gammaln(nmk) - nmk*log(2*pi);
iota = ones(n,1);
z = epe0*iota - 2*detval(:,1)*epe0d + detval(:,1).*detval(:,1)*eped;
% add quadratic terms based on prior for beta
Q1 = zeros(n,1);
Q2 = zeros(n,1);
xpxi = inv(xs'*xs);
sTI = sige*TI;
xpxis = inv(xs'*xs + sTI);
logdetx = log(det(xpxis));
C = C - 0.5*logdetx;
          for i=1:n;
           rho = detval(i,1);
           D = speye(nobs) - rho*W;
           bhat = xpxi*(xs'*D*ys);
           beta = xpxis*(xs'*D*ys + sTI*c); 
           Q1(i,1) = (c - beta)'*sTI*(c - beta);
           Q2(i,1) = (bhat - beta)'*(xs'*xs)*(bhat - beta);
          end;

den = C + detval(:,2) - nmk*log(z + Q1 + Q2);
out = real(den);

function [total,direct,indirect] = calc_sdm_effects(maxorder1,trs,strs,ntrs,rho,bhat,p)
% function to compute effects for 1 rho, 1 vector beta

%this forms the matrix of exponents needed to form G
saa=kron(ones(maxorder1,1),strs)';
saav=saa(:);
saavsub=saav(1:(ntrs*maxorder1));
ree1=reshape(saavsub,ntrs,maxorder1)';
op=ones(p,1);
trbig=trs(:,op)';

    arparmri=rho;
   
    %forming P
    blockbsdm=reshape(bhat(2:end,1),p,maxorder1);
       
    %forming G
    ree=arparmri.^ree1;   
    reblock=ree(1:maxorder1,1:maxorder1);
    ree(1:maxorder1,1:maxorder1)=reblock-tril(reblock)+(diag(ones(maxorder1,1)));
    
    %forming PG
    pg=blockbsdm*ree;%pg is also the total impacts by order (p by ntrs)
   
    %added contribution for iteration i
    pgtrbig=pg.*trbig;%direct
    pginbig=pg-pgtrbig;%indirect
    
total = pg;
direct = pgtrbig;
indirect = pginbig;

function [total,direct,indirect] = calc_sar_effects(maxorder1,trs,strs,ntrs,rho,bhat,p)
% function to compute effects for 1 rho, 1 vector beta

%this forms the matrix of exponents needed to form G
saa=kron(ones(maxorder1,1),strs)';
saav=saa(:);
saavsub=saav(1:(ntrs*maxorder1));
ree1=reshape(saavsub,ntrs,maxorder1)';
op=ones(p,1);
trbig=trs(:,op)';

    arparmri=rho;
   
    %forming P
    blockbsdm=reshape(bhat(2:end,1),p,maxorder1);
       
    %forming G
    ree=arparmri.^ree1;   
    reblock=ree(1:maxorder1,1:maxorder1);
    ree(1:maxorder1,1:maxorder1)=reblock-tril(reblock)+(diag(ones(maxorder1,1)));
    
    %forming PG
    pg=blockbsdm*ree;%pg is also the total impacts by order (p by ntrs)
   
    %added contribution for iteration i
    pgtrbig=pg.*trbig;%direct
    pginbig=pg-pgtrbig;%indirect
    
total = pg;
direct = pgtrbig;
indirect = pginbig;


function bounds = cr_interval(adraw,hperc)
% PURPOSE: Computes an hperc-percent credible interval for a vector of MCMC draws
% --------------------------------------------------------------------
% Usage: bounds = cr_interval(draws,hperc);
% where draws = an ndraw by nvar matrix
%       hperc = 0 to 1 value for hperc percentage point
% --------------------------------------------------------------------
% RETURNS:
%         bounds = a 1 x 2 vector with 
%         bounds(1,1) = 1-hperc percentage point
%         bounds(1,2) = hperc percentage point
%          e.g. if hperc = 0.95
%          bounds(1,1) = 0.05 point for 1st vector in the matrix
%          bounds(1,2) = 0.95 point  for 1st vector in the matrix
%          bounds(2,1) = 0.05 point for 2nd vector in the matrix
%          bounds(2,2) = 0.05 point for 2nd vector in the matrix
%          ...
% --------------------------------------------------------------------

% Written by J.P. LeSage

% This function takes a vector of MCMC draws and calculates
% an hperc-percent credible interval
[ndraw,ncols]=size(adraw);
botperc=round((0.50-hperc/2)*ndraw);
topperc=round((0.50+hperc/2)*ndraw);
bounds = zeros(ncols,2);
for i=1:ncols;
temp = sort(adraw(:,i),1);
bounds(i,:) =[temp(topperc,1) temp(botperc,1)];
end;


