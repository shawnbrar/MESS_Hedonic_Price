function results = sar_g(y,x,W,ndraw,nomit,prior)
% PURPOSE: Bayesian estimates of the spatial autoregressive model
%          y = rho*W*y + X*beta + e, e = N(0,sige*I_n), 
%          no priors for beta 
%          no priors for sige
%          no prior for rho
%-------------------------------------------------------------
% ATTENTION This functions modifies the sar_g function in the way the log-marginal
% likelihood is computed. Here, the log-marginal is computed in the same
% way as for the sar_conv_g function

% USAGE: results = sar_gnew(y,x,W,ndraw,nomit,prior)
% where: y = dependent variable vector (nobs x 1)
%        x = independent variables matrix (nobs x nvar), 
%            the intercept term (if present) must be in the first column of the matrix x
%        W = (nobs,nobs) weight matrix
%    ndraw = # of draws (use lots of draws, say 25,000 to 50,000
%    nomit = # of initial draws omitted for burn-in  (probably around 5,000
%    prior = a structure variable with:
%   
%            prior.ntraces = # of traces to use for Taylor's series (default=10)
%                            increasing this will lengthen the time
%                            required for trace calculations
%            prior.thin  = a thinning parameter for use in analyzing
%                          posterior distributions, default = 1 (no thinning of draws)
%                          recommended value for ndraw > 20,000 is 10
%                          default = 1
%            NOTE: thin is NOT used to determine how many times to MH-MC sample the
%            log-posterior using Monte Carlo integration, which is sampled
%            (ndraw-nomit) times
%            prior.sym   = 0,1 for trace calculations, 0 for non-symmetric
%            W-matrices, 1 for symmetric weight matrices (or similar to a
%            symmetric weight matrix), default = 0
%-------------------------------------------------------------
% RETURNS:  a structure:
%          results.meth     = 'sar_gnew'
%          results.beta     = posterior mean of bhat based on draws
%          results.rho      = posterior mean of rho based on draws
%          results.sige     = posterior mean of sige based on draws
%                             where m is the number of weight matrices used on input
%          results.nobs   = # of observations
%          results.nvar   = # of variables in x-matrix
%          results.p      = # of variables in x-matrix (excluding constant
%                           term if used)
%          results.beta_std = std deviation of beta draws
%          results.sige_std = std deviation of sige draws
%          results.rho_std  = std deviation of rho draws
%          results.sigma    = posterior mean of sige based on (e'*e)/(n-k)
%          results.bdraw    = bhat draws (1:thin:ndraw-nomit x nvar)
%          results.pdraw    = rho  draws (1:thin:ndraw-nomit x 1)
%          results.sdraw    = sige draws (1:thin:ndraw-nomit x 1)
%          results.thin     = thinning value from input
%          results.total    = a (1:thin:ndraw-nomit,p) total x-impacts
%          results.direct   = a (1:thin:ndraw-nomit,p) direct x-impacts
%          results.indirect = a (1:thin:ndraw-nomit,p) indirect x-impacts
%          results.lmarginal= a scalar log-marginal likelihood, from MH-MC
%                             (Metropolis-Hastings Monte Carlo) integration of the log-posterior
%          results.drawpost = a vector of log-posterior draws (ndraw-nomit)x1
%          results.rho_post = (ndraw-nomit) x 1 contains
%                              rhodraws used to evaluate the log-posterior
%          results.rho_mode   = modal value of rho from the log-posterior
%          results.logC_sar   = constants associated with log-marginal
%                               logC_sar = gammaln(dof) - dof*log(2*pi)  -0.5*lndetx_sar
%                               dof = (n - m)/2; lndetx_sar = log(det(xp'*x));
%          results.logm_profile = a profile of the log-margainal over [rho_post isum];
%                                 where: isum = exp(logp -adj) + adj; adj = max(logp)
%          results.ndraw  = # of draws
%          results.nomit  = # of initial draws omitted
%          results.y      = y-vector from input (nobs x 1)
%          results.yhat   = mean of posterior predicted (nobs x 1)
%          results.resid  = residuals, based on posterior means
%          results.rsqr   = r-squared based on posterior means
%          results.rbar   = adjusted r-squared%          results.ntraces = ntraces from input (or default = 10)
%          results.setup_time    = time for setup calculations 
%          results.sampling_time = time for MCMC sampling
%          results.effects_time  = time to calculate effects estimates
%          results.trace_time    = time to calculate traces for taylor approximations
%          results.time = sum of the above times
%          results.rmax   = 1  
%          results.rmin   = -1       
%          results.tflag  = 'plevel' (default) for printing p-levels
%                         = 'tstat' for printing bogus t-statistics 
%          results.cflag  = 1 for intercept term, 0 for no intercept term
%          results.tflag  = 'plevel' (default) for printing p-levels
%                         = 'tstat' for printing bogus t-statistics 
% --------------------------------------------------------------
% NOTES: - the intercept term (if you have one)
%          must be in the first column of the matrix x
% --------------------------------------------------------------
% SEE ALSO: (house_sar_demo.m, house_sar_demo2.m demos) 
% --------------------------------------------------------------
% REFERENCES: LeSage and Pace (2009) Chapter 5 on Bayesian estimation 
%             of spatial regression models.
% For lndet information see: Chapter 4 

%----------------------------------------------------------------

% written by:
% James P. LeSage, last updated 7/2020
% Dept of Finance & Economics
% Texas State University-San Marcos
% 601 University Drive
% San Marcos, TX 78666
% james.lesage@txstate.edu

% error checking on inputs
[n junk] = size(y);
[nc, k] = size(x);
if (nc ~= n)
       error('sar_g: wrong sized x-matrix');
end; 

results.nobs  = n;
results.nvar  = k;
results.novi=0; % Homoskedastic case
results.y = y; 
results.taylor=1;
    [n1,n2] = size(W);
    

if nargin == 5
    % use default arguments
    thin = 1;
    rmin = -1;     % use -1,1 rho interval as default
    rmax = 1;
    rho = 0.7; % starting values
    sige = 1;
    results.taylor = 1;
    results.rmin = -1;
    results.rmax = 1;
    results.thin=1;
    ccmin = 0.4;
    ccmax = 0.6;
    sym = 0;
    ntrs = 10;
    
    
elseif nargin == 6
     [rho,sige,rmin,rmax,ntrs,thin,ccmin,ccmax] = sar_parse(prior);

     results.thin = thin;
     results.rmin = rmin;
     results.rmax = rmax;
     results.ccmin = ccmin;
     results.ccmax = ccmax;
     

else
    error('sar_g: wrong # of input arguments to sar_g');
end

% check if the user handled the intercept term okay
    n = length(y);
    if sum(x(:,1)) ~= n
        tst = sum(x); % we may have no intercept term
        ind = find(tst == n); % we do have an intercept term
        if length(ind) > 0
            error('sar_g: intercept term must be in first column of the x-matrix');
        elseif length(ind) == 0 % case of no intercept term
            cflag = 1;
            pp = size(x,2);
        end;
    elseif sum(x(:,1)) == n % we have an intercept in the right place
        cflag = 0;
        pp = size(x,2)-1;
    end;
     
    results.cflag = cflag;
    results.p = pp;
    
    results.rmin = rmin;
    results.rmax = rmax;


% storage for draws
bsave = zeros(ndraw,k);
psave = zeros(ndraw,1);
ssave = zeros(ndraw,1);
acc_rate = zeros(ndraw,1);
cc_save = zeros(ndraw,1);
drawpost = zeros(ndraw-nomit,1);
rhopost = zeros(ndraw-nomit,1); % holds rho, beta, sige
sigpost = zeros(ndraw-nomit,1);

% ====== initializations
% compute this stuff once to save time

cc = 0.1;
acc = 0;
gacc = 0;

timet = clock; % start the timer
Wy = W*y;
xpx = x'*x;
Wty = [y Wy];
bd = xpx\(x'*Wty);

time = etime(clock,timet);
results.setup_time = time;

[traces,ntrs,ctime] = calc_trace_approx(W,ntrs);

results.trace_time = ctime;

timet = clock; % start the timer
    
for iter=1:ndraw
    
    % update beta
    AI = (xpx)\eye(k);
    
        tmp = [1
              -rho];
          
    Wys = Wty*tmp;

    b = x'*Wys;
    b0 = AI*b;
    bhat = norm_rnd(sige*AI) + b0;
    xb = x*bhat;
    
    % update sige
    V = sum((Wys -xb).^2)/2;
    sige=1/gamrand(n/2,V);

    
    % ====================== 
    % M-H sample rho
    % ====================== 
% Anonymous function
   lp_rho = @(r) cond_rho(r,Wty,bd,x,traces,ntrs,n);   % evaluates rho-conditional

    % obtain random-walk (tuned) proposed rho values
    rho2 = rho + cc*randn(1,1); % proposal for rho
    
    accept = 0;
    while accept == 0
        if ((rho2 > rmin) && (rho2 < rmax))
            accept = 1;
        else
            rho2 = rho + cc*randn(1,1);
        end
    end
   
    alpMH =  lp_rho(rho2) - lp_rho(rho);

    if alpMH > log(rand)
        rho = rho2;
        acc = acc + 1;
        
    end
        
    acc_rate(iter,1) = acc/iter;
    
    % update cc based on std of rho draws
    if acc_rate(iter,1) < ccmin
        cc = cc/1.1;
    end
    if acc_rate(iter,1) > ccmax
        cc = cc*1.1;
    end

    % save draws
    bsave(iter,:) = bhat';
    ssave(iter,1) = sige;
    psave(iter,1) = rho;
    
% Anonymous function    
    log_post = @(r) joint_post(r,Wty,bd,x,traces,ntrs,n); % evaluates log posterior for both rho and gamma

    if iter > nomit
            logpost = log_post(rho);
            drawpost(iter-nomit,1) = logpost;
            rhopost(iter-nomit,1) = rho;
            tmp = [1
                   -rho];
            btmp = bd*tmp;
            betapost(iter-nomit,:) = btmp;
            Wtys = Wty*tmp;
            xb = x*(bd*tmp);
            sigpost(iter-nomit,1) = ((Wtys - xb)'*(Wtys - xb))/n;
    end
    
%     tt=1:iter;
%     subplot(3,1,1),
%     plot(tt,bsave(1:iter,:));
%     subplot(3,1,2),
%     plot(tt,psave(1:iter,1));
%     subplot(3,1,3),
%     plot(tt,gsave(1:iter,:));
%     drawnow;
%     
    
    
end % end of draws loop
    
time = etime(clock,timet);
results.sampling_time = time;

results.cc = cc_save;
results.acc_rate = acc_rate;

results.thin = thin;
results.bdraw = bsave(nomit+1:thin:ndraw,:);
results.pdraw = psave(nomit+1:thin:ndraw,1);
results.sdraw = ssave(nomit+1:thin:ndraw,1);
results.acc_rate = acc_rate(nomit+1:ndraw,1);
results.drawpost = drawpost; % we don't want to thin these
results.rhopost = rhopost; % rho, beta, sige
results.betapost = betapost;
results.sigpost = sigpost;

% calculate log-marginal likelihood (using Mh-MC integration)
logp = results.drawpost;
post = [results.rhopost results.betapost results.sigpost];
[adj,mind] = max(logp); 
results.posterior_mode = post(mind,:);
isum = exp(logp -adj);
lndetx_sar = log(det(xpx));
% constant terms
dof = (n - 1)/2; 
D = (1 - 1/rmin); % from uniform prior on rho
logC_sar = -log(D) + gammaln(dof) - dof*log(2*pi)  -0.5*lndetx_sar;

results.logmarginal = mean(logp) + logC_sar;
results.logC_sar = logC_sar; % return constants
results.logm_profile = [post isum];


% compute posterior means for return arguments
bmean = mean(results.bdraw);
rho_mean = mean(results.pdraw);
smean = mean(results.sdraw);

results.sigma = smean;
results.beta = bmean';
results.rho = rho_mean;

% calculate fit statistics using posterior means for gamma
[nobs,nvar] = size(x);

        tmp = [1
            -rho_mean];
        
        %    ew = R*tmp;
        Wys = Wty*tmp;
        xb = x*bmean';
        
        ew = (Wys - xb);
        
        epe = ew'*ew;

ym = y - mean(y);
rsqr1 = epe;
rsqr2 = ym'*ym;
results.rsqr = 1- rsqr1/rsqr2; % r-squared
rsqr1 = rsqr1/(nobs-nvar);
rsqr2 = rsqr2/(nobs-1.0);
results.rbar = 1 - (rsqr1/rsqr2); % rbar-squared

% calculate effects estimates 
timet = clock; % start the timer

% =======================================
% here we use a trick for speed
% Use stochastic trace estimates based on
% mean(gamma), 
% then fill-in 2nd 3rd and 4th order traces
% that allow for variation in gamma based on MCMC draws

uiter=500;
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
    
pdraw = results.pdraw;
bdraw = results.bdraw;

niter = length(pdraw);

    total = zeros(niter,pp,101);
    direct = zeros(niter,pp,101);
    indirect = zeros(niter,pp,101);
    
    
    for iter=1:niter
        
        traces(1)=0;
        
        trs=[1;traces];
        ntrs=length(trs);
        trbig=trs';
        
        ree = 0:1:ntrs-1;
        
        %     rmat = zeros(1,ntrs);
        rmat = pdraw(iter,1).^ree;
        for j=1:pp
            if cflag == 0
                b = bdraw(iter,j+1);
            elseif cflag == 1
                b = bdraw(iter,j);
            end
            total(iter,j,:) = b*rmat;
            direct(iter,j,:) = (b*trbig).*rmat;
            indirect(iter,j,:) = total(iter,j,:) - direct(iter,j,:);
        end
        
    end
time = etime(clock,timet);
results.effects_time = time;
results.time = results.setup_time + results.effects_time + results.sampling_time + results.trace_time;
##
##% ====================================================================
##
total_out = zeros(niter,pp);
direct_out = zeros(niter,pp);

indirect_out = zeros(niter,pp);
for i=1:pp;
tmp = squeeze(total(:,i,:)); % an ndraw by 1 by ntraces matrix
total_out(:,i) = (sum(tmp'))'; % an ndraw by 1 vector
total_mean=mean(total_out)';
tmp = squeeze(indirect(:,i,:)); % an ndraw by 1 by ntraces matrix
indirect_out(:,i) = (sum(tmp'))'; % an ndraw by 1 vector
indirect_mean=mean(indirect_out)';
tmp = squeeze(direct(:,i,:)); % an ndraw by 1 by ntraces matrix
direct_out(:,i) = (sum(tmp'))'; % an ndraw by 1 vector
direct_mean=mean(direct_out)';
end;

results.total_mean = total_mean;
results.direct_mean = direct_mean;
results.indirect_mean = indirect_mean;
    
results.total = total_out;
results.direct = direct_out;
results.indirect = indirect_out;



results.meth  = 'sar_g';
results.ndraw = ndraw;
results.nomit = nomit;
results.tflag = 'plevel';


% =========================================================================
% support functions below
% =========================================================================

function [logp] = joint_post(rho,Wty,bd,x,traces,ntrs,n)
% PURPOSE: evaluate the  joint distribution of rho and gamma
%  spatial autoregressive model using sparse matrix algorithms
% ---------------------------------------------------
%  USAGE: cout = joint_post(rho,Wy,H,T2,T3,T4,n)
%  where:  rho  = spatial autoregressive parameter
%          Wy   = AR dependence vector [y W1*y W2*y .. WL*y]
%          HWy  = Hat matrix times Wy matrix
%         traces = traces used to find det(I-rho*W) 
%                 using Taylor series approximation 
%          n    =  nobs
%    
% ---------------------------------------------------
%  RETURNS: a joint distribution value used in Monte Carlo integration
%           to produce the log-marginal likelihood
%  --------------------------------------------------

tmp = [1
     -rho];
 
Wtys = Wty*tmp;
xb = x*(bd*tmp);
ew = (Wtys - xb);

epe = ew'*ew;


    
        traces(1)=0;
        
        ree = 1:1:ntrs;
        
        tt=1:ntrs;
        rmat = rho.^ree./tt;
        
        tmp = rmat*traces;
    
    lndet = -tmp;


logp =  lndet - (n/2)*log(epe);



function rhoc = cond_rho(rho,Wty,bd,x,traces,ntrs,n)
% PURPOSE: evaluate the  conditional distribution of rho 
%  spatial autoregressive model using sparse matrix algorithms
% ---------------------------------------------------
%  USAGE:cout = cond_rho(rho,Wy,xb,sige,aTa,aTTa,aTTTa,n)
%  where:  rho  = spatial autoregressive parameter
%          Wy   = AR dependence vector [y W1*y W2*y .. WL*y]
%          HWy  = Hat matrix times Wy matrix
%         traces=  traces used to find det(I-rho*W) 
%                 using Taylor series approximation 
%          n    =  nobs
%       
% ---------------------------------------------------
%  RETURNS: a conditional used in Metropolis-Hastings sampling
%   NOTE: Used only for sar_conv_g
%  --------------------------------------------------

tmp = [1
     -rho];

 Wtys = Wty*tmp;
xb = x*(bd*tmp);
ew = (Wtys - xb);

epe = ew'*ew;
                   
    ree = 1:1:ntrs;
    
    tt=1:ntrs;
    
    rmat = (rho.^ree)./tt;
    
    tmp = rmat*traces;

    lndet = -tmp;
    


rhoc =  lndet -(n/2)*log(epe);


% ===========================================================================


function [rho,sige,rmin,rmax,ntrs,thin,ccmin,ccmax] = sar_parse(prior)

% PURPOSE: parses input arguments for sar_conv_g models
% ---------------------------------------------------
%  USAGE: [rval,rho,sige,rmin,rmax,novi_flag] =  sar_parse(prior,k)
% returns values set by user or default values 
% ---------------------------------------------------

% written by:
% James P. LeSage, last updated 1/2018
% Dept of Finance & Economics
% Texas State University-San Marcos
% 601 University Drive
% San Marcos, TX 78666
% james.lesage@txstate.edu


% set defaults
rmin = -1;     % use -1,1 rho interval as default
rmax = 1;
rho = 0.7; % starting values
sige = 1;
thin = 1; % default to no thinning
ccmin = 0.4;
ccmax = 0.6;
ntrs = 10;
fields = fieldnames(prior);
nf = length(fields);
if nf > 0
    for i=1:nf
        if strcmp(fields{i},'thin')
            thin = prior.thin;
        elseif strcmp(fields{i},'ccmin')
            ccmin = prior.ccmin;
        elseif strcmp(fields{i},'ccmax')
            ccmax = prior.ccmax;
        elseif strcmp(fields{i},'ntrs')
            ntrs = prior.ntrs;
        end
    end
    
    
else % the user has input a blank info structure
    % so we use the defaults
end


function [traces,ntrs,ctime] = calc_trace_approx(W,ntrs)
% PURPOSE: calculate ntrs order trace matrices for Taylor Series
%  approximation to the log-determinant

timet = clock;
[n,nm] = size(W);


maxorderu=ntrs;
nobs = n;
rv=randn(nobs,1);
tracew=zeros(maxorderu,1);
wjjju=rv;
wpow = W;
for jjj=1:maxorderu
    tmp = wjjju'*wpow*wjjju;
    tracew(jjj)=mean(mean(tmp));
    wpow = wpow*W;
end

traces=tracew;

ntrs=length(traces);
traces(1)=0;



ctime = etime(clock,timet);


function x=gamrand(alpha,lambda)
% Gamma(alpha,lambda) generator using Marsaglia and Tsang method
% Algorithm 4.33
if alpha>1
    d=alpha-1/3; c=1/sqrt(9*d); flag=1;
    while flag
        Z=randn;
        if Z>-1/c
            V=(1+c*Z)^3; U=rand;
            flag=log(U)>(0.5*Z^2+d-d*V+d*log(V));
        end
    end
    x=d*V/lambda;
else
    x=gamrand(alpha+1,lambda);
    x=x*rand^(1/alpha);
end



