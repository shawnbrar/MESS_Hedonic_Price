function prt_sar_conv_bma_g(results,vnames,fid,fmt)
% PURPOSE: Prints output using SAR results structures
%---------------------------------------------------
% USAGE: prt_sar_conv(results,vnames,fid)
% Where: results = a structure returned by a sar_conv_g model
%        vnames  = an optional vector of variable names
%        fid     = optional file-id for printing results to a file
%                  (defaults to the MATLAB command window)
%        fmt     = format string, e.g., '%12.4f' (default)
%--------------------------------------------------- 
%  NOTES: e.g. vnames = strvcat('y','const','x1','x2');
%         e.g. fid = fopen('ols.out','wr');
%  use prt_spat(results,[],fid) to print to a file with no vnames               
% --------------------------------------------------
%  RETURNS: nothing, just prints the sar_conv_g estimation results
% --------------------------------------------------
% SEE ALSO: prt, plt
%---------------------------------------------------   

% written by:
% James P. LeSage, 4/2018
% Dept of Finance & Economics
% Texas State University-San Marcos
% 601 University Drive
% San Marcos, TX 78666
% jlesage@spatial-econometrics.com

if ~isstruct(results)
 error('bma_prt_sar_conv requires structure argument');
elseif nargin == 1
 nflag = 0; fid = 1; fmt = '%12.4f';
elseif nargin == 2
 fid = 1; nflag = 1; fmt = '%12.4f';
elseif nargin == 3
 nflag = 0;  fmt = '%12.4f';
 [vsize junk] = size(vnames); % user may supply a blank argument
   if vsize > 0
   nflag = 1;          
   end;
elseif nargin == 4
    fmt = fmt
else
 error('Wrong # of arguments to bma_prt_sar_conv');
end;


nvar = results.nvar;
nobs = results.nobs;
cflag = results.cflag;
nmat = results.nmat; % Number of considered connectivy matrices.

if (nflag == 1) % the user supplied variable names
[tst_n nsize] = size(vnames);
 if tst_n ~= nvar+1
 fprintf(fid,'Wrong # of variable names in bma_prt_sar_conv -- check vnames argument \n');
 fprintf(fid,'will use generic variable names \n');
 nflag = 0;
 end
end;

% handling of vnames
Vname = 'Variable';
if nflag == 0 % no user-supplied vnames or an incorrect vnames argument
    if cflag == 1 % a constant term

        Vname = strvcat(Vname,'constant');
     for i=1:nvar-1
        tmp = ['variable ',num2str(i)];
        Vname = strvcat(Vname,tmp);
     end;
 
    elseif cflag == 0 % no constant term

     for i=1:nvar
        tmp = ['variable ',num2str(i)];
        Vname = strvcat(Vname,tmp);
     end;
    end;
 
     
% add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
    
% add gamma parameter names
for ii=1:nmat
        Vname = strvcat(Vname,['gamma' num2str(ii)]);
end;


elseif (nflag == 1) % the user supplied variable names
    if cflag == 1 % no constant term
    Vname = 'Variable';
     for i=1:nvar
        Vname = strvcat(Vname,vnames(i+1,:));
     end;
    % add spatial rho parameter name
        Vname = strvcat(Vname,'rho');
        % add gamma parameter names
for ii=1:nmat
        Vname = strvcat(Vname,['gamma' num2str(ii)]);
end;

     elseif cflag == 0 % a constant term
     Vname = 'Variable';
     for i=1:nvar
        Vname = strvcat(Vname,vnames(i+1,:));
     end;
    % add spatial rho parameter name
        Vname = strvcat(Vname,'rho');
        % add gamma parameter names
for ii=1:nmat
        Vname = strvcat(Vname,['gamma' num2str(ii)]);
end;

    end; % end of cflag issue       
 
end; % end of nflag issue


% extract posterior means
bout = [results.beta
        results.rho
        results.gamma];
    
% sige = results.sige;
    tmp1 = std(results.bdraw);
    tmp2 = std(results.pdraw);
    tmp3 = std(results.gdraw);
    
    bstd = [tmp1'
            tmp2
            tmp3']; 
        
tout = bout./bstd;

if strcmp(results.tflag,'tstat')
    tstat = bout./bstd;
    % find t-stat marginal probabilities
    tout = tdis_prb(tstat,results.nobs);
    results.tstat = bout./bstd; % trick for printing below
else % find plevels
    bout = plims(results.bdraw);
    pout = plims(results.pdraw);
    gout = plims(results.gdraw);
    bout = [bout
            pout
            gout];


end

% rsqr = results.rsqr;

% do effects estimates
% =======================================================
% a set of draws for the effects/impacts distribution
     direct_out =  plims(results.direct);
     indirect_out =  plims(results.indirect);
     total_out =  plims(results.total);
        

fprintf(fid,'\n');
fprintf(fid,'Bayesian Model Average of spatial autoregressive convex W models \n');
if (nflag == 1)
fprintf(fid,'Dependent Variable  = %16s \n',vnames(1,:));
end;
fprintf(fid,'BMA Log-marginal    = %9.4f \n',results.lmarginal);
% cstats2 = chainstats(results.drawpost);
% fprintf(fid,'Log-marginal MCerror= %9.6f\n',cstats2(1,2));
% fprintf(fid,'R-squared           = %9.4f \n',rsqr);
% fprintf(fid,'Rbar-squared        = %9.4f \n',results.rbar);
% fprintf(fid,'mean of sige draws  = %9.4f \n',mean(results.sdraw));
% fprintf(fid,'posterior mode sige = %9.4f \n',results.sig_mode);
fprintf(fid,'Nobs, Nvars         = %6d,%6d \n',results.nobs,results.nvar);
fprintf(fid,'# weight matrices   = %6d \n', results.nmat);
fprintf(fid,'ndraws,nomit        = %6d,%6d \n',results.ndraw,results.nomit);
% fprintf(fid,'total time in secs  = %9.4f   \n',results.time);
fprintf(fid,'total time          = %9.4f \n',results.time);
% fprintf(fid,'time for sampling   = %9.4f \n',results.sampling_time);
% fprintf(fid,'time for Taylor     = %9.4f \n',results.taylor_time);
fprintf(fid,'thinning for draws  = %6d   \n',results.thin);
fprintf(fid,'min and max rho     = %9.4f,%9.4f \n',results.rmin,results.rmax);

fprintf(fid,'***************************************************************\n');
nd = (results.ndraw-results.nomit)/results.thin;
fprintf(fid,'      MCMC diagnostics ndraws = %d \n',nd);
cstats=chainstats([results.bdraw results.pdraw results.gdraw]);
in.cnames = strvcat('Mean','MC error','tau','Geweke');
in.rnames = Vname;
in.fmt = strvcat('%12.4f','12.8f','%10.6f','%10.6f');
in.width = 10000;
% pmode = [results.beta_mode results.rho_mode results.gamma_mode]';
out = [cstats];
mprint(out,in);

fprintf(fid,'***************************************************************\n');
fprintf(fid,'      Posterior Estimates \n');

if strcmp(results.tflag,'tstat')
% now print coefficient estimates, t-statistics and probabilities
tout = norm_prb(results.tstat); % find asymptotic z (normal) probabilities
      
tmp = [bout bstd results.tstat tout];  % matrix to be printed
% column labels for printing results
bstring = 'Coefficient';sstring = 'Std dev'; tstring = 'Asymptot t-stat'; pstring = 't-probability';
cnames = strvcat(bstring,sstring,tstring,pstring);
in.cnames = cnames;
in.rnames = Vname;
in.fmt = fmt;
in.fid = fid;

mprint(tmp,in);
else % use p-levels for Bayesian results
% % column labels for printing results
cnames = strvcat('lower 0.01','lower 0.05','median','upper 0.95','upper 0.99');
in.cnames = cnames;
in.rnames = Vname;

in.fmt = fmt;
in.fid = fid;
mprint(bout,in);
end;

% now print x-effects estimates

if strcmp(results.tflag,'tstat')
bstring = 'Coefficient'; 
tstring = 't-stat'; 
pstring = 't-prob';
lstring = 'lower 01';
ustring = 'upper 99';
cnames = strvcat(bstring,tstring,pstring,lstring,ustring);
else
cnames = strvcat('lower 0.01','lower 0.05','median','upper 0.95','upper 0.99');
end;

ini.width = 2000;
ini.cnames = cnames;

% print effects estimates

if cflag == 0 % we have an intercept
vnameso = strvcat(Vname(3:end-1-nmat,:));
elseif cflag == 1
vnameso = strvcat(Vname(2:end-1-nmat,:));
end    

ini.rnames = strvcat('Direct',vnameso);
ini.fmt = fmt;
ini.fid = fid;

% set up print out matrix
printout = direct_out;
mprint(printout,ini);

printout = indirect_out;
ini.rnames = strvcat('Indirect',vnameso);
mprint(printout,ini);

printout = total_out;
ini.rnames = strvcat('Total',vnameso);
mprint(printout,ini);




function yp=plims(x,p)
%PLIMS Empirical quantiles
% plims(x,p)  calculates p quantiles from columns of x

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 1.5 $  $Date: 2012/09/27 11:47:39 $

if nargin<2
p = [0.005, 0.025, 0.5, 0.975, 0.995];   
%    p = [0.25,0.5,0.75];
end
[n,m] = size(x); if n==1; n=m;end
y = interp1(sort(x),(n-1)*p+1);
if m > 1
yp = y';
else
yp = y;
end

function stats=chainstats(chain)
%CHAINSTATS some statistics from the MCMC chain
% chainstats(chain,results)
%    chain    nsimu*npar MCMC chain
%    results  output from mcmcrun function

% $Revision: 1.4 $  $Date: 2009/08/13 15:47:35 $

mcerr = bmstd(chain)./sqrt(size(chain,1));

[z,p]  = geweke(chain);
tau    = iact(chain);
stats  = [mean(chain)',mcerr',tau', p'];



function [z,p]=geweke(chain,a,b)
%GEWEKE Geweke's MCMC convergence diagnostic
% [z,p] = geweke(chain,a,b)
% Test for equality of the means of the first a% (default 10%) and
% last b% (50%) of a Markov chain.
% See:
% Stephen P. Brooks and Gareth O. Roberts.
% Assessing convergence of Markov chain Monte Carlo algorithms.
% Statistics and Computing, 8:319--335, 1998.

% ML, 2002
% $Revision: 1.3 $  $Date: 2003/05/07 12:22:19 $

[nsimu,npar]=size(chain);

if nargin<3
  a = 0.1;
  b = 0.5;
end

na = floor(a*nsimu);
nb = nsimu-floor(b*nsimu)+1;

if (na+nb)/nsimu >= 1
  error('Error with na and nb');
end

m1 = mean(chain(1:na,:));
m2 = mean(chain(nb:end,:));

%%% Spectral estimates for variance
sa = spectrum0(chain(1:na,:));
sb = spectrum0(chain(nb:end,:));

z = (m1-m2)./(sqrt(sa/na+sb/(nsimu-nb+1)));
p = 2*(1-nordf(abs(z)));


function s=bmstd(x,b)
%BMSTD standard deviation calculated from batch means
% s = bmstd(x,b) - x matrix - b length of the batch
% bmstd(x) gives an estimate of the Monte Carlo std of the 
% estimates calculated from x

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 1.2 $  $Date: 2012/09/27 11:47:34 $

[n,p] = size(x);

if nargin<2
  b = max(10,fix(n/20));
end

inds = 1:b:(n+1);
nb = length(inds)-1;
if nb < 2
  error('too few batches');
end

y = zeros(nb,p);

for i=1:nb
  y(i,:)=mean(x(inds(i):inds(i+1)-1,:));
end

% calculate the estimated std of MC estimate
s = sqrt( sum((y - repmat(mean(x),nb,1)).^2)/(nb-1)*b );

function s=spectrum0(x)
%SPECTRUM0 Spectral density at frequency zero
% spectrum0(x) spectral density at zero for columns of x

% ML, 2002
% $Revision: 1.3 $  $Date: 2003/05/07 12:22:19 $

[m,n]= size(x);
s = zeros(1,n);
for i=1:n
  spec = spectrum(x(:,i),m);
  s(i) = spec(1);
end


function [y,f]=spectrum(x,nfft,nw)
%SPECTRUM Power spectral density using Hanning window
%  [y,f]=spectrum(x,nfft,nw) 

% See also: psd.m in Signal Processing Toolbox 

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 1.4 $  $Date: 2012/09/27 11:47:40 $

if nargin < 2 | isempty(nfft)
  nfft = min(length(x),256);
end
if nargin < 3 | isempty(nw)
  nw = fix(nfft/4);
end
noverlap = fix(nw/2);

% Hanning window
w = .5*(1 - cos(2*pi*(1:nw)'/(nw+1)));
% Daniel
%w = [0.5;ones(nw-2,1);0.5];
n = length(x);
if n < nw
    x(nw)=0;  n=nw;
end
x = x(:);

k = fix((n-noverlap)/(nw-noverlap)); % no of windows
index = 1:nw;
kmu = k*norm(w)^2; % Normalizing scale factor
y = zeros(nfft,1);
for i=1:k
% xw = w.*detrend(x(index),'linear');
  xw = w.*x(index);
  index = index + (nw - noverlap);
  Xx = abs(fft(xw,nfft)).^2;
  y = y + Xx;
end

y = y*(1/kmu); % normalize

n2 = floor(nfft/2);
y  = y(1:n2);
f  = 1./n*(0:(n2-1));


function y=nordf(x,mu,sigma2)
% NORDF the standard normal (Gaussian) cumulative distribution.
% NORPF(x,mu,sigma2) x quantile, mu mean, sigma2 variance

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 1.6 $  $Date: 2012/09/27 11:47:38 $

if nargin < 2, mu     = 0; end
if nargin < 3, sigma2 = 1; end

%y = 0.5*erf(-Inf,sqrt(2)*0.5*x);
y = 0.5+0.5*erf((x-mu)/sqrt(sigma2)/sqrt(2));


function [tau,m] = iact(dati)
%IACT estimates the integrated autocorrelation time
%   using Sokal's adaptive truncated periodogram estimator.

% Originally contributed by Antonietta Mira by name sokal.m

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 1.2 $  $Date: 2012/09/27 11:47:37 $

if length(dati) == prod(size(dati))
  dati = dati(:);
end

[mx,nx] = size(dati);
tau = zeros(1,nx);
m   = zeros(1,nx);

x  = fft(dati);
xr = real(x);
xi = imag(x);
xr = xr.^2+xi.^2; %%%controllare questo
xr(1,:)=0;
xr=real(fft(xr));
var=xr(1,:)./length(dati)/(length(dati)-1);

for j = 1:nx
  if var(j) == 0
    continue
  end
  xr(:,j)=xr(:,j)./xr(1,j);
  sum=-1/3;
  for i=1:length(dati)
    sum=sum+xr(i,j)-1/6;
    if sum<0
      tau(j)=2*(sum+(i-1)/6);
      m(j)=i;
      break
    end
  end
end
