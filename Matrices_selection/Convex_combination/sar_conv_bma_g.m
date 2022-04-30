function results = sar_conv_bma_g(y,x,Wmatrices,ndraw,nomit,prior)
% PURPOSE: Bayesian model averaged estimates for all M combinations
%          of 2 or more W-matrices, spatial autoregressive models
%          using a convex combination of m different W-matrices
%          M = 2^L - L - 1, where L = # of weight matrices 
%          y = rho*Wc*y + X*beta + e, e = N(0,sige*I_n), 
%          Wc = g1*W1 + g2*W2 + ... + (1-g1-g2- ... -gL)*WL
%          no priors for beta 
%          no priors for sige
%          uniform (-1,1) prior for rho
%          uniform (0,1) prior for g1,g2, ... gm
% ATTENTION : If your Matlab licence does not allow for parallel computing,
% replace "parfor" on line 184 by "for"
%-------------------------------------------------------------
% USAGE: results = sar_conv_gbma(y,x,Wmatrices,ndraw,nomit,prior)
% where: y = dependent variable vector (nobs x 1)
%        x = independent variables matrix (nobs x nvar), 
%            the intercept term (if present) must be in the first column of the matrix x
%        Wmatrices = (nobs,L*nobs)
%        e.g., Wmatrices = [W1 W2 ... WL]
%        where each W1, W2, ... WL are (nobs x nobs) row-normalized weight matrices
%    ndraw = # of draws (use lots of draws, say 25,000 to 50,000
%    nomit = # of initial draws omitted for burn-in  (probably around 5,000
%    prior = a structure variable with:
%            prior.torder  = 4 for 4th order Taylor series, 5 for 5th order
%            prior.thin  = a thinning parameter for use in analyzing
%                          posterior distributions, default = 1 (no thinning of draws)
%                          recommended value for ndraw > 20,000 is 10
%                          default = 1
%            NOTE: thin is NOT used to determine how many times to MH-MC sample the
%            log-posterior using Monte Carlo integration, which is sampled
%            (ndraw-nomit) times
%-------------------------------------------------------------
% PRINTS out M log-marginals and M model probabilities, and
% RETURNS a structure variable: 
%          a structure:
%          results.meth     = 'sar_conv_gbma'
%          results.beta     = posterior mean of BMA bhat based on draws
%          results.rho      = posterior mean of BMA rho based on draws
%          results.sige     = posterior mean of BMA sige based on draws
%          results.gamma    = L x 1 vector of posterior means for  BMA g1,g2, ... gL
%                             where L is the number of weight matrices used on input
%          results.gamma_all = an M x L matrix of gamma estimates for all M models
%          results.nobs   = # of observations
%          results.nvar   = # of variables in x-matrix
%          results.p      = # of variables in x-matrix (excluding constant term if used)
%          results.nmat   = # of W-matrices from input
%          results.sigma    = posterior mean of BMA sige based on (e'*e)/(n-k)
%          results.bdraw    = bhat BMA draws (1:thin:ndraw-nomit x nvar)
%          results.pdraw    = rho  BMA draws (1:thin:ndraw-nomit x 1)
%          results.sdraw    = sige BMA draws (1:thin:ndraw-nomit x 1)
%          results.gdraw    = gamma BMA draws (1:thin:ndraw-nomit x L)
%          results.thin     = thinning value from input
%          results.total    = a (1:thin:ndraw-nomit,p) BMA total x-impacts
%          results.direct   = a (1:thin:ndraw-nomit,p) BMA direct x-impacts
%          results.indirect = a (1:thin:ndraw-nomit,p) BMA indirect x-impacts
%          results.lmarginal= a scalar log-marginal BMA likelihood, from MH-MC
%                             (Metropolis-Hastings Monte Carlo) integration of the log-posterior
%          results.prob     = Mx1 vector of model probabilities
%          results.ndraw  = # of draws
%          results.nomit  = # of initial draws omitted
%          results.y      = y-vector from input (nobs x 1)
%          results.time   = time for MCMC sampling
%          results.rmax   =  1  
%          results.rmin   = -1       
%          results.cflag  = 0 for intercept term, 1 for no intercept term
% --------------------------------------------------------------
% NOTES: - the intercept term (if you have one)
%          must be in the first column of the matrix x
% --------------------------------------------------------------
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
       error('sar_conv_bma_g: wrong sized x-matrix');
end

if ndraw <= 1000
       error('sar_conv_bma_g: ndraw<=1000, increase ndraw to at least 10000');
end


results.nobs  = n;
results.nvar  = k;
results.y = y; 
% check if the user handled the intercept term okay
    n = length(y);
    if sum(x(:,1)) ~= n
        tst = sum(x); % we may have no intercept term
        ind = find(tst == n); % we do have an intercept term
        if length(ind) > 0
            error('sar_conv_bma_g: intercept term must be in first column of the x-matrix');
        elseif length(ind) == 0 % case of no intercept term
            cflag = 1;
            pp = size(x,2);
        end
    elseif sum(x(:,1)) == n % we have an intercept in the right place
        cflag = 0;
        pp = size(x,2)-1;
    end
    prior1=prior;
    prior1.pp=pp;
     prior1.cflag=cflag;
    results.cflag = cflag;
    results.p = pp;
    prior1.plt=0;

    [n1,n2] = size(Wmatrices);
    nmat = n2/n1;
    if nmat~= round(nmat)
        error('sar_conv_bma_g: wrong sized W-matrices');
    elseif n1 ~= n
        error('sar_conv_bma_g: wrong sized W-matrices');
    elseif nmat < 2
        error('sar_conv_bma_g: only one W-matrix');
    end
    
    
results.meth = 'sar_conv_bma_g';
results.nmat=nmat; 
v = 1:nmat;
cnt = 0;
for i=1:nmat
    tmp=nchoosek(v,i);
    [n1,~] = size(tmp);
 cnt = cnt + n1;   
    
end


v = 1:nmat;
cnt = 0;
for i=1:nmat
    tmp=nchoosek(v,i);
    [n1,n2] = size(tmp);
 cnt = cnt + n1;   
    
end;


out = zeros(cnt,nmat);
cnt = 0;
for i=1:nmat
    tmp=nchoosek(v,i);
    [n1,n2] = size(tmp);
    out(1+cnt:n1+cnt,1:n2) = tmp;
 cnt = cnt + n1;   
    
end



% ==========================================
indices = out(nmat+1:end,:);

np = length(indices);

% form a structure variable with the sets of matrices
% outside the parfor loop
for i = 1:np
    ind = find(indices(i,:) > 0);
    [~,nw] = size(ind);
    wi = indices(i,ind);

    Wtmp = [];

    for j=1:nw
         begi=(wi(j)-1)*n+1;
         endi=wi(j)*n;
        Wtmp = [Wtmp Wmatrices(:,begi:endi)];
    end
    Wmatrix(i).model = Wtmp;
    
end

timet=clock;
 parfor iter = 1:np
    
    res = sar_conv_g(y,x,Wmatrix(iter).model,ndraw,nomit,prior1);
    
    temp(iter) = res;
 end
 time=etime(clock,timet);
 results.time=time;
thin=temp(1).thin;
results.rmin=temp(1).rmin;
results.rmax=temp(1).rmax;
 if thin ~= 1
    neff = (ndraw-nomit)/prior.thin;
 else
     neff = ndraw-nomit;
 end
 
bsave = zeros(np,neff,k);

psave = zeros(np,neff);
ssave = zeros(np,neff);
smean=zeros(np,1);
gsave = zeros(np,neff,nmat);
gprint = zeros(np,nmat);
dsave = zeros(np,neff,pp);
isave = zeros(np,neff,pp);
tsave = zeros(np,neff,pp);
rho_mean = zeros(np,1);
beta_mean = zeros(np,k);
% g_mean=zeros(np,nmat);
logm = zeros(np,1);
ttime = zeros(np,1);

for ii=1:np
    
    
    logm(ii,1) = temp(ii).logmarginal;
    rho_mean(ii,1) = temp(ii).rho;
    beta_mean(ii,:)=temp(ii).beta';
%     g_mean(ii,:)=temp(ii).gamma';
    dsave(ii,:,:) = temp(ii).direct;
    isave(ii,:,:) = temp(ii).indirect;
    tsave(ii,:,:)=temp(ii).total;
    bsave(ii,:,:) = temp(ii).bdraw;
    psave(ii,:) = temp(ii).pdraw';
    ssave(ii,:) = temp(ii).sdraw';
    smean(ii,1)=temp(ii).sige;
    ttime(ii,1)=temp(ii).time;
    position = [];
    for j=1:nmat
        pos = find(indices(ii,:)==j);
        if length(pos) > 0
            position = [position j];
        else
            position = [position 0];
        end
    end
    cnt = 1;
    for j=1:nmat
        if position(1,j) > 0
            gsave(ii,:,position(1,j)) = temp(ii).gdraw(:,cnt);
            gprint(ii,position(1,j)) = temp(ii).gamma(cnt,1);
            cnt = cnt+1;
        end
    end
    
end 

probs = model_probs(logm);
 prt_out3 = [logm probs rho_mean gprint];

cnames = strvcat('logm','Prob','rho');
for i=1:nmat
    cnames = strvcat(cnames,['W' num2str(i)]);
end

in.cnames = cnames;

% rnames = strvcat('Models');
rnames = strvcat('Models');
for i=1:np
    rnames = strvcat(rnames,['Model ' num2str(i)]); 
end
% for i = 1:nmat
%     rnames=strvcat(rnames,['W' num2str(i)]);
% end

in.rnames = rnames;

in.width = 10000;
fmt = strvcat('%12.3f');
% for i=1:nmat;
%     fmt = strvcat(fmt,'%5d');
% end;
in.fmt = fmt;


hi=find(probs==max(probs));

bma_rho=probs'*rho_mean;
bma_sige = probs'*smean;
bma_gamma = (probs'*gprint);
bma_logm = sum(probs.*logm);
prt_out3 = [prt_out3
            bma_logm  1 bma_rho bma_gamma 
            logm(hi,1) probs(hi,:) rho_mean(hi,:) gprint(hi,:)];
 in.rnames = strvcat(rnames, 'BMA','highest');
mprint(prt_out3,in);

%% Computation of the BMA values for coefficients and effects 

    bma_beta=(probs'*beta_mean)';

%    bma_gamma=sum(matmul(probs,gamma_mean))'; 

results.probs=probs;
results.rho=bma_rho;
results.beta=bma_beta;
results.gamma=bma_gamma';
results.sige = bma_sige;


bma_b=zeros(neff,k);
bma_dir=zeros(neff,pp);
bma_ind=zeros(neff,pp);
bma_tot=zeros(neff,pp);
bma_r=zeros(neff,1);
bma_s=zeros(neff,1);
bma_g=zeros(neff,nmat);

for j=1:k
    bma_b(:,j)=sum(matmul(probs,bsave(:,:,j)));
end
    bma_r(:,1)=sum(matmul(probs,psave));
    bma_s(:,1)=sum(matmul(probs,ssave));
for j = 1:nmat
    bma_g(:,j)=sum(matmul(probs,gsave(:,:,j)));
end
for j = 1:pp
    bma_dir(:,j)=sum(matmul(probs,dsave(:,:,j)));
    bma_ind(:,j)=sum(matmul(probs,isave(:,:,j)));
    bma_tot(:,j)=sum(matmul(probs,tsave(:,:,j)));
end

results.gamma_all=gprint;
results.total=bma_tot;
results.direct=bma_dir;
results.indirect=bma_ind;
results.meth= 'sar_conv_bma_g';
results.tflag = 'plevel';
results.ndraw=ndraw;
results.nomit=nomit;

results.bdraw=bma_b;
results.pdraw=bma_r;
results.sdraw=bma_s;
results.gdraw=bma_g;
results.lmarginal=bma_logm;
results.thin=thin;







