function results = sar_conv_g(y,x,Wmatrices,ndraw,nomit,prior)
% PURPOSE: Bayesian estimates of the spatial autoregressive model
%          using a convex combination of m different W-matrices
%          y = rho*Wc*y + X*beta + e, e = N(0,sige*I_n), 
%          Wc = g1*W1 + g2*W2 + ... + (1-g1-g2- ... -gm)*Wm
%          no priors for beta 
%          no priors for sige
%          uniform (-1,1) prior for rho
%          uniform [0,1] prior for g1,g2, ... gm
%-------------------------------------------------------------
% USAGE: results = sar_conv_g(y,x,Wmatrices,ndraw,nomit,prior)
% where: y = dependent variable vector (nobs x 1)
%        x = independent variables matrix (nobs x nvar), 
%            the intercept term (if present) must be in the first column of the matrix x
%        Wmatrices = (nobs,m*nobs)
%        e.g., Wmatrices = [W1 W2 ... Wm]
%        where each W1, W2, ... Wm are (nobs x nobs) row-normalized weight matrices
%    ndraw = # of draws (use lots of draws, say 25,000 to 50,000
%    nomit = # of initial draws omitted for burn-in  (probably around 5,000
%    prior.plt = 1 for plotting of MCMC draws, 0 for no plots, default = 0
%    prior = a structure variable with:
%            prior.thin  = a thinning parameter for use in analyzing
%                          posterior distributions, default = 1 (no thinning of draws)
%                          recommended value for ndraw > 20,000 is 10
%                          default = 1
%            NOTE: thin is NOT used to determine how many times to MH-MC sample the
%            log-posterior using Monte Carlo integration, which is sampled
%            (ndraw-nomit) times
%    prior.T2 = pre-calculated 2nd order Taylor series traces fed to the function
%    prior.T3 = pre-calculated 2nd order Taylor series traces fed to the function
%    prior.T4 = pre-calculated 2nd order Taylor series traces fed to the function
% these will be calculated by this function, but can be pre-calculated and
% fed to the function in the case of Monte Carlo studies to save some time
%-------------------------------------------------------------
% RETURNS:  a structure:
%          results.meth     = 'sar_conv_g'
%          results.beta     = posterior mean of bhat based on draws
%          results.rho      = posterior mean of rho based on draws
%          results.sige     = posterior mean of sige based on draws
%          results.gamma    = m x 1 vector of posterior means for g1,g2, ... gm
%                             where m is the number of weight matrices used on input
%          results.nobs   = # of observations
%          results.nvar   = # of variables in x-matrix
%          results.p      = # of variables in x-matrix (excluding constant
%                           term if used)
%          results.beta_std = std deviation of beta draws
%          results.sige_std = std deviation of sige draws
%          results.rho_std  = std deviation of rho draws
%          results.g_std    = m x 1 vector of posterior std deviations
%          results.sigma    = posterior mean of sige based on (e'*e)/(n-k)
%          results.bdraw    = bhat draws (1:thin:ndraw-nomit x nvar)
%          results.pdraw    = rho  draws (1:thin:ndraw-nomit x 1)
%          results.sdraw    = sige draws (1:thin:ndraw-nomit x 1)
%          results.gdraw    = gamma draws (1:thin:ndraw-nomit x m)
%          results.thin     = thinning value from input
%          results.total    = a (1:thin:ndraw-nomit,p) total x-impacts
%          results.direct   = a (1:thin:ndraw-nomit,p) direct x-impacts
%          results.indirect = a (1:thin:ndraw-nomit,p) indirect x-impacts
%          results.lmarginal= a scalar log-marginal likelihood, from MH-MC
%                             (Metropolis-Hastings Monte Carlo) integration of the log-posterior
%          results.drawpost = a vector of log-posterior draws (ndraw-nomit)x1
%          results.betapost = a matrix of log-posterior draws for beta (ndraw-nomit x k)
%          results.sigpost  = a vector of log-posterior draws for sige (ndraw-nomit x 1)
%          results.rho_gamma = (ndraw-nomit) x m+1 contains
%                              [rhodraws gammdraws] used to evaluate the log-posterior
%          results.rho_mode   = modal value of rho from the log-posterior
%          results.gamma_mode = modal values of gamma vector from the log-posterior
%          results.beta_mode  = modal values of beta vector from the log-posterior
%          results.sig_mode   = modal values of beta vector from the log-posterior
%          results.logC_sar   = constants associated with log-marginal
%                               logC_sar = gammaln(dof) - dof*log(2*pi)  -0.5*lndetx_sar
%                               dof = (n - m)/2; lndetx_sar = log(det(xp'*x));
%          results.logm_profile = a profile of the log-margainal over [rho_gamma isum];
%                                 where: isum = exp(logp -adj); adj = max(logp)
%          results.ndraw  = # of draws
%          results.nomit  = # of initial draws omitted
%          results.y      = y-vector from input (nobs x 1)
%          results.yhat   = mean of posterior predicted (nobs x 1)
%          results.resid  = residuals, based on posterior means
%          results.rsqr   = r-squared based on posterior means
%          results.rbar   = adjusted r-squared
%          results.taylor = taylor from input;
%          results.sampling_time = time for MCMC sampling
%          results.effects_time  = time to calculate effects estimates
%          results.trace_time    = time to calculate traces for taylor/chebyshev approximations
%          results.rmax   =  1  
%          results.rmin   = -1       
%          results.cflag  = 0 for intercept term, 1 for no intercept term
%          results.T2 = 2nd order Taylor series traces
%          results.T3 = 3rd order Taylor series traces
%          results.T4 = 4th order Taylor series traces
% --------------------------------------------------------------
% NOTES: - the intercept term (if you have one)
%          must be in the first column of the matrix x
% --------------------------------------------------------------
% SEE ALSO: (sar_conv_g_demo.m demos) 
% --------------------------------------------------------------
% REFERENCES: Debarsy and LeSage (2020) 
% Bayesian model averaging for spatial autoregressive models
% based on convex combinations of different types of connectivity
% matrices, Journal of Business & Economics Statistics,
% Forthcoming.


% error checking on inputs
[n junk] = size(y);
[nc, k] = size(x);
if (nc ~= n)
       error('sar_conv_g: wrong sized x-matrix');
end

if ndraw <= 1000
       error('sar_conv_g: ndraw<=1000, increase ndraw to at least 10000');
end


results.nobs  = n;
results.nvar  = k;
results.y = y; 

    [n1,n2] = size(Wmatrices);
    m = n2/n1;
    if m~= round(m)
        error('sar_conv_g: wrong sized W-matrices');
    elseif n1 ~= n
        error('sar_conv_g: wrong sized W-matrices');
    elseif m < 2
        error('sar_conv_g: only one W-matrix');
    end;
   results.nmat=m;
if nargin == 5
    % use default arguments
    thin = 1;
    rmin = -1;     % use -1,1 rho interval as default
    rmax = 1;
    rho = 0.7; % starting values
    sige = 1;
    gamma = ones(m,1)*(1/m);
    results.rmin = -1;
    results.rmax = 1;
    ccmin = 0.4;
    ccmax = 0.6;
    ddmin = 0.1;
    ddmax = 0.4;
    tr_flag=0;
    
     results.ccmin = ccmin;
     results.ccmax = ccmax;
     results.ddmin = ddmin;
     results.ddmax = ddmax;
    
elseif nargin == 6

     [rho,sige,rmin,rmax,gamma,thin,ccmin,ccmax,ddmin,ddmax,T2,T3,T4,tr_flag,plt_flag] = sar_parse(prior,m);

% check the thinning parameter
 eff = (ndraw -nomit)/thin;
if eff < 999
    warning('sar_conv_g: < 1000 draws after thining');
end

     thin = round(thin);
     results.thin = thin;
     results.rmin = rmin;
     results.rmax = rmax;
     results.ccmin = ccmin;
     results.ccmax = ccmax;
     results.ddmin = ddmin;
     results.ddmax = ddmax;
     

else
    error('sar_conv_g: wrong # of input arguments to sar_conv_g');
end

% check if the user handled the intercept term okay
    n = length(y);
    if sum(x(:,1)) ~= n
        tst = sum(x); % we may have no intercept term
        ind = find(tst == n); % we do have an intercept term
        if length(ind) > 0
            error('sar_conv_g: intercept term must be in first column of the x-matrix');
        elseif length(ind) == 0 % case of no intercept term
            cflag = 1;
            pp = size(x,2);
        end
    elseif sum(x(:,1)) == n % we have an intercept in the right place
        cflag = 0;
        pp = size(x,2)-1;
    end
    results.cflag = cflag;
    results.p = pp;
      
results.rmin = rmin;
results.rmax = rmax;

% storage for draws
bsave = zeros(ndraw,k);
psave = zeros(ndraw,1);
ssave = zeros(ndraw,1);
gsave = zeros(ndraw,m);
drawpost = zeros(ndraw-nomit,1);
rho_gamma = zeros(ndraw-nomit,m+1);
betapost = zeros(ndraw-nomit,k);
sigpost = zeros(ndraw-nomit,1);

   % ====== initializations
ccmin = 0.4;
ccmax = 0.6;
ddmin = 0.10;
ddmax = 0.40;
acc_rate = zeros(ndraw,1);
gacc_rate = zeros(ndraw,1);
cc_save = zeros(ndraw,1);
gflag = zeros(ndraw,1);
cc = 0.1;
acc = 0;
gacc = 0;
dd = 3.0;
dd_save = zeros(ndraw,1);

Wy = y;
xpx = x'*x;

begi = 1;
endi = n;
for ii=1:m
    Wy = [Wy Wmatrices(:,begi:endi)*y];
    begi = begi + n;
    endi = endi + n;
end

bd = xpx\(x'*Wy);
    
AI = (xpx)\eye(k);

if (tr_flag == 0)
[T2,T3,T4,ctime] = calc_taylor_approx4(Wmatrices);
results.taylor_time = ctime;
elseif (tr_flag == 1)
results.taylor_time = 0;
end


results.T2=T2;
results.T3=T3;
results.T4=T4;

timet = clock; % start the timer


    noo = 1000;
    
for iter=1:ndraw
    
    % update beta
    
        tmp = [1
              -rho*gamma];
          
    Wys = Wy*tmp;

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
   lp_rho = @(r) cond_rho4(r,gamma,Wy,bd,x,T2,T3,T4,n,k);   % evaluates rho-conditional on gamma
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
    
    ru = unif_rnd(1,0,1);
    
    if alpMH > 0
        p = 1;
    else
        ratio = exp(alpMH);
        p = min(1,ratio);
    end
    if (ru < p)
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

    if cc < 0.001
    cc = 0.1;
    end
    if cc > 1.000
    cc = 0.1;
    end
    
    cc_save(iter,1) = cc;
    % ====================== 
    % M-H sample gamma
    % ====================== 
% Anonymous function
    lp_gamma = @(g) cond_gamma4(g,rho,Wy,bd,x,T2,T3,T4,n,k); % evaluates gamma conditional on rho
    if iter < noo
        % obtain a block of gamma candidates using reversible jump
        gtst = zeros(m-1,1);
        % flip a coin
        coin = unif_rnd(m-1,0,1);
        for jj=1:m-1
            if coin(jj,1) <= (1/3)
                gtst(jj,1) = unif_rnd(1,0,gamma(jj,1));
            elseif (coin(jj,1)  > (1/3) && coin(jj,1)  <= (2/3))
                gtst(jj,1) = gamma(jj,1);
            elseif coin(jj,1)  > (2/3)
                gtst(jj,1) = unif_rnd(1,gamma(jj,1),1);
            end
        end
        
        gnew = [gtst
            1-sum(gtst)];
        accept = 0;
        cntr = 1;
        while (accept == 0)
            if all(gnew >= 0)
                accept = 1;
            else
                gtst = zeros(m-1,1);
                % flip a 3-headed coin
                coin = unif_rnd(m-1,0,1);
                for jj=1:m-1
                    if coin(jj,1) <= (1/3)
                        gtst(jj,1) = unif_rnd(1,0,gamma(jj,1));
                    elseif (coin(jj,1)  > (1/3) && coin(jj,1)  <= (2/3))
                        gtst(jj,1) = gamma(jj,1);
                    elseif coin(jj,1)  > (2/3)
                        gtst(jj,1) = unif_rnd(1,gamma(jj,1),1);
                    end
                end
                gnew = [gtst
                    1-sum(gtst)];
                cntr = cntr+1;
            end
        end
        
    elseif iter == noo
        rmp = std(gsave(iter-1000+1:iter-1,:));
        gam_std = rmp';
        ind = find(gam_std < 0.001);
        if length(ind) > 0
            gam_std(ind,1) = 0.01;
        end
                
        % obtain a block of gamma candidates using reversible jump
        gtst = zeros(m-1,1);
        % flip a coin
        coin = unif_rnd(m-1,0,1);
        for jj=1:m-1
            if coin(jj,1) <= (1/3)
                gtst(jj,1) = unif_rnd(1,gamma(jj,1) - dd*gam_std(jj,1),gamma(jj,1));
            elseif (coin(jj,1)  > (1/3) && coin(jj,1)  <= (2/3))
                gtst(jj,1) = gamma(jj,1);
            elseif coin(jj,1)  > (2/3)
                gtst(jj,1) = unif_rnd(1,gamma(jj,1),gamma(jj,1) + dd*gam_std(jj,1));
            end
        end
        
        gnew = [gtst
            1-sum(gtst)];
        accept = 0;
        cntr = 1;
        while (accept == 0)
            if all(gnew >= 0)
                accept = 1;
            else
                gtst = zeros(m-1,1);
                % flip a 3-headed coin
                coin = unif_rnd(m-1,0,1);
                for jj=1:m-1
                    if coin(jj,1) <= (1/3)
                        gtst(jj,1) = unif_rnd(1,gamma(jj,1) - dd*gam_std(jj,1),gamma(jj,1));
                    elseif (coin(jj,1)  > (1/3) && coin(jj,1)  <= (2/3))
                        gtst(jj,1) = gamma(jj,1);
                    elseif coin(jj,1)  > (2/3)
                        gtst(jj,1) = unif_rnd(1,gamma(jj,1),gamma(jj,1) + dd*gam_std(jj,1));
                    end
                end
                gnew = [gtst
                    1-sum(gtst)];
                cntr = cntr+1;
            end
        end

    else % use gamma std to produce draws
                
        % obtain a block of gamma candidates using reversible jump
        gtst = zeros(m-1,1);
        % flip a coin
        coin = unif_rnd(m-1,0,1);
        for jj=1:m-1
            if coin(jj,1) <= (1/3)
                gtst(jj,1) = unif_rnd(1,gamma(jj,1) - dd*gam_std(jj,1),gamma(jj,1));
            elseif (coin(jj,1)  > (1/3) && coin(jj,1)  <= (2/3))
                gtst(jj,1) = gamma(jj,1);
            elseif coin(jj,1)  > (2/3)
                gtst(jj,1) = unif_rnd(1,gamma(jj,1),gamma(jj,1) + dd*gam_std(jj,1));
            end
        end
        
        gnew = [gtst
            1-sum(gtst)];
        accept = 0;
        cntr = 1;
        while (accept == 0)
            if all(gnew >= 0)
                accept = 1;
            else
                gtst = zeros(m-1,1);
                % flip a 3-headed coin
                coin = unif_rnd(m-1,0,1);
                for jj=1:m-1
                    if coin(jj,1) <= (1/3)
                        gtst(jj,1) = unif_rnd(1,gamma(jj,1) - dd*gam_std(jj,1),gamma(jj,1));
                    elseif (coin(jj,1)  > (1/3) && coin(jj,1)  <= (2/3))
                        gtst(jj,1) = gamma(jj,1);
                    elseif coin(jj,1)  > (2/3)
                        gtst(jj,1) = unif_rnd(1,gamma(jj,1),gamma(jj,1) + dd*gam_std(jj,1));
                    end
                end
                gnew = [gtst
                    1-sum(gtst)];
                cntr = cntr+1;
            end
        end
    end
    
    % evaluate conditional for the block of gamma proposals
       
        alpMH =  lp_gamma(gnew) - lp_gamma(gamma);
        
        ru = unif_rnd(1,0,1);
    
        if alpMH > 0
            p = 1;
        else
            ratio = exp(alpMH);
            p = min(1,ratio);
        end
    
    if (ru < p)
        gflag(iter,1) = 1;
        gacc = gacc + 1;
        gamma = gnew;
    end
        
        gacc_rate(iter,1) = gacc/iter;
        
        % update dd based on std of gamma draws
        if gacc_rate(iter,1) < ddmin
            dd = dd/1.1;
        end
        if gacc_rate(iter,1) > ddmax
            dd = dd*1.1;
        end
        
        if dd>3.0
            dd=3.0;
        elseif dd<1.0
            dd=1.0;
        end
        
        
     dd_save(iter,1) = dd;
     
    % save draws
    bsave(iter,:) = bhat';
    ssave(iter,1) = sige;
    psave(iter,1) = rho;
    gsave(iter,:) = gamma';
    
% Anonymous function    
    log_post = @(r,g) joint_post4(r,g,Wy,bd,x,T2,T3,T4,n,k); % evaluates log posterior for both rho and gamma


    if ( mod(iter, noo) ==0 )
        rmp = std(gsave(iter-noo+1:iter,:));
        gam_std = rmp';
        ind = find(gam_std < 0.001);
        if length(ind) > 0
            gam_std(ind,1) = 0.01;
        end
     
        if plt_flag == 1
        subplot(2,2,1),       
        tt=iter-noo+1:iter-1;
        subplot(2,2,1),
        plot(tt,gsave(iter-noo+1:iter-1,:));
        xlabel('gammas draws');
        subplot(2,2,2),
        plot(tt,psave(iter-noo+1:iter-1,1));
        xlabel('rho draws');
        subplot(2,2,3),
        plot(tt,bsave(iter-noo+1:iter-1,:));
        xlabel('beta draws');
        subplot(2,2,4),
        plot(tt,ssave(iter-noo+1:iter-1,:));
        xlabel('sigma draws');
        drawnow;
        end
        
    end
    
    
    if iter > nomit
        
            logpost = log_post(rho,gamma);
            drawpost(iter-nomit,1) = logpost;
            rho_gamma(iter-nomit,:) = [rho gamma'];
            tmp = [1
                   -rho*gamma];
            btmp = bd*tmp;
            betapost(iter-nomit,:) = btmp';
            sigpost(iter-nomit,1) = ((Wy*tmp - x*btmp)'*(Wy*tmp - x*btmp))/n;
    end
    
%     
    
    
end % end of draws loop
    
time = etime(clock,timet);
results.sampling_time = time;

results.gflag = gflag;

results.gacc_rate = gacc_rate;
results.cc = cc_save;
results.dd = dd_save;
results.acc_rate = acc_rate;

results.thin = thin;
results.bdraw = bsave(nomit+1:thin:ndraw,:);
results.pdraw = psave(nomit+1:thin:ndraw,1);
results.sdraw = ssave(nomit+1:thin:ndraw,1);
results.gdraw = gsave(nomit+1:thin:ndraw,:);
% results.acc_rate = acc_rate(nomit+1:ndraw,1);
% results.gacc_rate = gacc_rate(nomit+1:ndraw,1);
results.drawpost = drawpost; % we don't want to thin these
results.rho_gamma = rho_gamma;
results.betapost = betapost;
results.sigpost = sigpost;

% calculate log-marginal likelihood (using Mh-MC integration)
logp = results.drawpost;
rho_gamma = results.rho_gamma;
[adj,mind] = max(logp); 
results.rho_mode = rho_gamma(mind,1);
results.beta_mode = betapost(mind,:);
results.sig_mode = sigpost(mind,1);
results.gamma_mode = rho_gamma(mind,2:end);
isum = exp(logp -adj);
lndetx_sar = log(det(xpx));
% constant terms

dof = (n - m)/2; % we must include the # of weight matrices
D = (1 - 1/rmin); % from uniform prior on rho
logC_sar = -log(D) + gammaln(dof) - dof*log(2*pi)  -0.5*lndetx_sar;

results.logmarginal = mean(logp) + logC_sar;
results.logC_sar = logC_sar; % return constants
results.logm_profile = [rho_gamma betapost sigpost isum];
% results.logm_profile = [rho_gamma isum];


% compute posterior means for return arguments
bmean = mean(results.bdraw);
rho_mean = mean(results.pdraw);
smean = mean(results.sdraw);
gmean = mean(results.gdraw);

results.sige = smean;
results.beta = bmean';
results.rho = rho_mean;
results.gamma = gmean';

% calculate fit statistics using posterior means for gamma


[nobs,nvar] = size(x);

tmp=[1
    -rho_mean*gmean'];
Wys=Wy*tmp;

        
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

uiter=150;
maxorderu=100;
nobs = n;

begi = 1;
endi = n;
Wch = sparse(n,n);
for ii=1:m;
    Wch = Wch + gmean(1,ii)*sparse(Wmatrices(:,begi:endi));
    begi = begi + n;
    endi = endi + n;
end;

Wch = sparse(Wch);

    rv=randn(nobs,uiter);
    tracew=zeros(maxorderu,1);
    wjjju=rv;
    for jjj=1:maxorderu
        wjjju=Wch*wjjju;
        tracew(jjj)=mean(mean(rv.*wjjju));
        
    end
    
    traces=[tracew];
    
pdraw = results.pdraw;
gdraw = results.gdraw;
bdraw = results.bdraw;

niter = length(pdraw);

    total = zeros(niter,pp,101);
    direct = zeros(niter,pp,101);
    indirect = zeros(niter,pp,101);
    
    for iter=1:niter;
        
        traces(1)=0;
        ghat = gdraw(iter,:)';
        a2=kron(ghat,ghat);
        traces(2) = a2'*T2/nobs;
        a3=kron(a2,ghat);
        traces(3) = a3'*T3/nobs;
        a4=kron(a3,ghat);
        traces(4) = a4'*T4/nobs;
        
        trs=[1;traces];
        ntrs=length(trs);
        trbig=trs';
        
        ree = 0:1:ntrs-1;
        
        %     rmat = zeros(1,ntrs);
        rmat = pdraw(iter,1).^ree;
        for j=1:pp;
            if cflag == 1
                b = bdraw(iter,j);
            elseif cflag == 0
                b = bdraw(iter,j+1);
            end;
            total(iter,j,:) = b*rmat;
            direct(iter,j,:) = (b*trbig).*rmat;
            indirect(iter,j,:) = total(iter,j,:) - direct(iter,j,:);
        end;
        
    end;
time = etime(clock,timet);
results.effects_time = time;

results.time = results.effects_time + results.sampling_time + results.taylor_time;

% ====================================================================

total_out = zeros(niter,pp);
direct_out = zeros(niter,pp);
indirect_out = zeros(niter,pp);
for i=1:pp;
tmp = squeeze(total(:,i,:)); % an ndraw by 1 by ntraces matrix
total_out(:,i) = (sum(tmp'))'; % an ndraw by 1 vector
tmp = squeeze(indirect(:,i,:)); % an ndraw by 1 by ntraces matrix
indirect_out(:,i) = (sum(tmp'))'; % an ndraw by 1 vector
tmp = squeeze(direct(:,i,:)); % an ndraw by 1 by ntraces matrix
direct_out(:,i) = (sum(tmp'))'; % an ndraw by 1 vector
end;

results.total = total_out;
results.direct = direct_out;
results.indirect = indirect_out;


results.meth  = 'sar_conv_g';
results.ndraw = ndraw;
results.nomit = nomit;
results.tflag = 'plevel';


% =========================================================================
% support functions below
% =========================================================================

function [logp] = joint_post4(rho,gamma,Wy,bd,x,T2,T3,T4,n,k)
% PURPOSE: evaluate the  joint distribution of rho and gamma
%  spatial autoregressive model using sparse matrix algorithms
% ---------------------------------------------------
%  USAGE: cout = joint_post(rho,gamma,Wy,bd,x,T2,T3,T4,n)
%  where:  rho  = spatial autoregressive parameter
%        gamma  = convex combination parameters
%          Wy   = AR dependence vector [y W1*y W2*y .. WL*y]
%          bd  = xpx\(x'*Wy);
%         T2, T3, T4 = traces used to find det(I-rho*W) 
%          n    =  nobs
% ---------------------------------------------------
%  RETURNS: a joint distribution value used in Monte Carlo integration
%           to produce the log-marginal likelihood
%  --------------------------------------------------

tmp = [1
     -rho*gamma];
Wys = Wy*tmp;
xb = x*(bd*tmp);
ew = (Wys - xb);

epe = ew'*ew;


    a2=kron(gamma,gamma);
    aTa = a2'*T2;
    
    a3=kron(a2,gamma);
    aTTa =a3'*T3;
    
    a4=kron(a3,gamma);
    aTTTa =a4'*T4;
    
    wAw4 = (rho*rho*aTa/2) + (rho*rho*rho*aTTa/3) + (rho*rho*rho*rho*aTTTa/4);
    
    lndet = -wAw4;
    

logp =  lndet - ((n-k)/2)*log(epe);



function rhoc = cond_rho4(rho,gamma,Wy,bd,x,T2,T3,T4,n,k)
% PURPOSE: evaluate the  conditional distribution of rho given gamma
%  spatial autoregressive model using sparse matrix algorithms
% ---------------------------------------------------
%  USAGE:cout = cond_rho(rho,gamma,Wy,bd,x,T2,T3,T4,n)
%  where:  rho  = spatial autoregressive parameter
%        gamma  = convex combination parameters
%          Wy   = AR dependence vector [y W1*y W2*y .. WL*y]
%          bd  = xpx\(x'*Wy);
%         T2, T3, T4 = traces used to find det(I-rho*W) 
%                 using Chebyshev or Taylor series approximation 
%          n    =  nobs
% ---------------------------------------------------
%  RETURNS: a conditional used in Metropolis-Hastings sampling
%  NOTE: called only by sar_conv_g
%  --------------------------------------------------

tmp = [1
     -rho*gamma];
 
Wys = Wy*tmp;
xb = x*(bd*tmp);
ew = (Wys - xb);

epe = ew'*ew;


    a2=kron(gamma,gamma);    
    aTa = a2'*T2;
    
    a3=kron(a2,gamma);
    aTTa =a3'*T3;
    
    a4=kron(a3,gamma);
    aTTTa =a4'*T4;
    
    wAw4 = (rho*rho*aTa/2) + (rho*rho*rho*aTTa/3) + (rho*rho*rho*rho*aTTTa/4);
    
    lndet = -wAw4;
    


rhoc =  lndet -((n-k)/2)*log(epe);



% ====================================================
function gamc = cond_gamma4(gamma,rho,Wy,bd,x,T2,T3,T4,n,k)
% PURPOSE: evaluate the  conditional distribution of rho given gamma
%  spatial autoregressive model using sparse matrix algorithms
% ---------------------------------------------------
%  USAGE:cout = cond_gamma(gamma,rho,Wy,,bd,x,T2,T3,T4,n)
%  where:  rho  = spatial autoregressive parameter
%        gamma  = convex combination parameters
%          Wy   = AR dependence vector [y W1*y W2*y .. WL*y]
%          bd  = xpx\(x'*Wy);
%         T2, T3, T4 = traces used to find det(I-rho*W) 
%          n    =  nobs
% ---------------------------------------------------
%  RETURNS: a conditional used in Metropolis-Hastings sampling
%  NOTE: called only by sar_conv_g
%  --------------------------------------------------

tmp = [1
     -rho*gamma];

Wys = Wy*tmp;
xb = x*(bd*tmp);
ew = (Wys - xb);

epe = ew'*ew;

    a2=kron(gamma,gamma);
    aTa = a2'*T2;
    
    a3=kron(a2,gamma);
    aTTa =a3'*T3;
    
    a4=kron(a3,gamma);
    aTTTa =a4'*T4;
    
    wAw4 = (rho*rho*aTa/2) + (rho*rho*rho*aTTa/3) + (rho*rho*rho*rho*aTTTa/4);
    
    lndet = -wAw4;
    

gamc =  lndet -((n-k)/2)*log(epe);

% ====================================================

% ===========================================================================


function [rho,sige,rmin,rmax,gamma,thin,ccmin,ccmax,ddmin,ddmax,T2,T3,T4,tr_flag,plt_flag] = sar_parse(prior,m)
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
gamma = ones(m,1)*(1/m);
rmin = -1;     % use -1,1 rho interval as default
rmax = 1;
rho = 0.7; % starting values
sige = 1;
thin = 1; % default to no thinning
ccmin = 0.4;
ccmax = 0.6;
ddmin = 0.1;
ddmax = 0.4;
tr_flag = 0;
plt_flag = 0;
T2 = [];
T3 = [];
T4 = [];


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
         elseif strcmp(fields{i},'ddmin')
            ddmin = prior.ddmin;
         elseif strcmp(fields{i},'ddmax')
            ddmax = prior.ddmax;          
        elseif strcmp(fields{i},'T2') && strcmp(fields{i},'T3') && strcmp(fields{i},'T4')
            T2 = prior.T2;    
            T3 = prior.T3;    
            T4 = prior.T4;    
            tr_flag = 1;
        elseif strcmp(fields{i},'plt')
            plt_flag = prior.plt;          
        end
    end
    
    
else % the user has input a blank info structure
    % so we use the defaults
end

function [tmat2,tmat3,tmat4,ctime] = calc_taylor_approx4(Wmatrices)
% PURPOSE: calculate 4th other trace matrices for Taylor Series
% approximation to the log-determinant

[n,nm] = size(Wmatrices);
m = nm/n;

% ========================================
% 2nd order
tmat2 = zeros(m*m,1);
% ========================================
% 3rd order
tmat3 = zeros(m*m*m,1);
% ========================================
% 4th order
tmat4 = zeros(m*m*m*m,1);


    tic;

    begi = 1;
    endi = n;
    cnti = 1;
    cntj = 1;
    cntk = 1;
    for ii=1:m;
        begj = 1;
        endj = n;
        Wi = sparse(Wmatrices(:,begi:endi));
        for jj=1:m;
            begk = 1;
            endk = n;
            Wj = sparse(Wmatrices(:,begj:endj));
            ijsave = (Wi*Wj);
            if (cnti <= m*m)
                tmat2(cnti,1) = sum(sum(Wi.*Wj'));
                cnti = cnti + 1;
            end;
            
            for kk=1:m;
                begl = 1;
                endl = n;
                Wk = sparse(Wmatrices(:,begk:endk));
                ijksave = ijsave*Wk;
                if (cntj <= m*m*m)
                    tmat3(cntj,1) = sum(sum((ijsave).*Wk'));
                    cntj = cntj + 1;
                end;
                for ll=1:m;
                    Wl = sparse(Wmatrices(:,begl:endl));
                    tmat4(cntk,1) = sum(sum((ijksave).*Wl'));
                    cntk = cntk+1;
                    begl = begl+n;
                    endl = endl+n;
                end;
                begk = begk + n;
                endk = endk + n;
            end;
            begj = begj+n;
            endj = endj+n;
        end;
        begi = begi + n;
        endi = endi + n;
    end;
    ctime = toc;
    


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




