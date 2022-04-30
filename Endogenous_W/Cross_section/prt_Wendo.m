function prt_Wendo(results,vnames,fid)
% PURPOSE: Prints results structures returned by most functions
%          by calling the appropriate printing function
%---------------------------------------------------
% USAGE: prt(results,vnames,fid)
% Where: results = a results structure returned an econometric function
%        vnames  = an optional vector of variable names
%        fid     = file-id for printing results to a file
%                  (defaults to the MATLAB command window)
%---------------------------------------------------               
%                 e.g. vnames = ['y    ',
%                                'x1   ',  NOTE: fixed width
%                                'x2   ',        like all MATLAB
%                                'cterm'];       strings
%                 e.g. fid = fopen('ols.out','wr');
% --------------------------------------------------
% NOTES: you may use prt(results,[],fid) to print
%        output to a file with no vnames
%        this is simply a wrapper function that calls another function
% --------------------------------------------------        
% RETURNS:
%        nothing, just prints the regression results
% --------------------------------------------------
% SEE ALSO: plt()
%---------------------------------------------------   

% written by:
% James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jlesage@spatial-econometrics.com

% error checking on inputs
if ~isstruct(results)
error('prt: requires a structure input');
elseif nargin == 3
arg = 0;
 [vsize junk] = size(vnames); % user may supply a blank argument
   if vsize > 0
   arg = 3;          
   end;
elseif nargin == 2
arg = 2;
elseif nargin == 1
arg = 1;
else
error('Wrong # of inputs to prt');
end;

method = results(1).meth;

% call appropriate printing routine
switch method


case {'sar_Wendo'}
     % call prt_sar
     if arg == 1
     prt_sar_Wendo(results);
     elseif arg == 2
     prt_sar_Wendo(results,vnames);
     elseif arg == 3
     prt_sar_Wendo(results,vnames,fid);
     else
     prt_sar_Wendo(results,[],fid);
     end;

 
    
otherwise
error('results structure not known by prt function');

end;






function prt_sar_Wendo(results,vnames,fid)
% PURPOSE: Prints output using SAR results structures
%---------------------------------------------------
% USAGE: prt_sar(results,vnames,fid)
% Where: results = a structure returned by a SAR model
%        vnames  = an optional vector of variable names
%        fid     = optional file-id for printing results to a file
%                  (defaults to the MATLAB command window)
%--------------------------------------------------- 
%  NOTES: e.g. vnames = strvcat('y','const','x1','x2');
%         e.g. fid = fopen('ols.out','wr');
%  use prt_spat(results,[],fid) to print to a file with no vnames               
% --------------------------------------------------
%  RETURNS: nothing, just prints the SAR results
% --------------------------------------------------
% SEE ALSO: prt, plt
%---------------------------------------------------   

% written by:
% James P. LeSage, 3/2010
% Dept of Finance & Economics
% Texas State University-San Marcos
% 601 University Drive
% San Marcos, TX 78666
% jlesage@spatial-econometrics.com

if ~isstruct(results)
 error('prt_sar requires structure argument');
elseif nargin == 1
 nflag = 0; fid = 1;
elseif nargin == 2
 fid = 1; nflag = 1;
elseif nargin == 3
 nflag = 0;
 [vsize junk] = size(vnames); % user may supply a blank argument
   if vsize > 0
   nflag = 1;          
   end;
else
 error('Wrong # of arguments to prt_sar');
end;


nvar = results.nvar;
nobs = results.nobs;
cflag = results.cflag;

if (nflag == 1) % the user supplied variable names
[tst_n nsize] = size(vnames);
 if tst_n ~= nvar+1
 fprintf(fid,'Wrong # of variable names in prt_sdm -- check vnames argument \n');
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

elseif (nflag == 1) % the user supplied variable names
    if cflag == 0 % no constant term
    Vname = 'Variable';
     for i=1:nvar
        Vname = strvcat(Vname,vnames(i+1,:));
     end;
    % add spatial rho parameter name
        Vname = strvcat(Vname,'rho');
     elseif cflag == 1 % a constant term
     Vname = 'Variable';
     for i=1:nvar
        Vname = strvcat(Vname,vnames(i+1,:));
     end;
    % add spatial rho parameter name
        Vname = strvcat(Vname,'rho');
    end; % end of cflag issue       
 
end; % end of nflag issue



switch results.meth

case {'sar_Wendo'} % <=================== max lik spatial autoregressive model

nobs = results.nobs;
nvar = results.nvar;
ndraw = results.ndraw;

% do effects estimates
% =======================================================
% a set of draws for the effects/impacts distribution
 total    = results.total;
 indirect = results.indirect;
 direct   = results.direct;

% Compute means, std deviation and upper and lower 0.99 intervals
iter = ndraw;
p = results.p;
total_out = zeros(p,5);
total_save = zeros(ndraw,p);
for i=1:p;
tmp = squeeze(total(:,i,:)); % an ndraw by 1 by ntraces matrix
total_mean = mean(tmp);
total_std = std(tmp);
% Bayesian 0.99 credible intervals
% for the cumulative total effects
total_sum = (sum(tmp'))'; % an ndraw by 1 vector
cum_mean = cumsum(mean(tmp));
cum_std = cumsum(std(tmp));
total_save(:,i) = total_sum;
bounds = cr_interval(total_sum,0.99);
cmean = mean(total_sum);
smean = std(total_sum);
ubounds = bounds(1,1);
lbounds = bounds(1,2);
total_out(i,:) = [cmean cmean./smean tdis_prb(cmean./smean,nobs) lbounds ubounds];
end;

% now do indirect effects
indirect_out = zeros(p,5);
indirect_save = zeros(ndraw,p);
for i=1:p;
tmp = squeeze(indirect(:,i,:)); % an ndraw by 1 by ntraces matrix
indirect_mean = mean(tmp);
indirect_std = std(tmp);
% Bayesian 0.95 credible intervals
% for the cumulative indirect effects
indirect_sum = (sum(tmp'))'; % an ndraw by 1 vector
cum_mean = cumsum(mean(tmp));
cum_std = cumsum(std(tmp));
indirect_save(:,i) = indirect_sum;
bounds = cr_interval(indirect_sum,0.99);
cmean = mean(indirect_sum);
smean = std(indirect_sum);
ubounds = bounds(1,1);
lbounds = bounds(1,2);
indirect_out(i,:) = [cmean cmean./smean tdis_prb(cmean./smean,nobs) lbounds ubounds  ];
end;


% now do direct effects
direct_out = zeros(p,5);
direct_save = zeros(ndraw,p);
for i=1:p;
tmp = squeeze(direct(:,i,:)); % an ndraw by 1 by ntraces matrix
direct_mean = mean(tmp);
direct_std = std(tmp);
% Bayesian 0.95 credible intervals
% for the cumulative direct effects
direct_sum = (sum(tmp'))'; % an ndraw by 1 vector
cum_mean = cumsum(mean(tmp));
cum_std = cumsum(std(tmp));
direct_save(:,i) = direct_sum;
bounds = cr_interval(direct_sum,0.99);
cmean = mean(direct_sum);
smean = std(direct_sum);
ubounds = bounds(1,1);
lbounds = bounds(1,2);
direct_out(i,:) = [cmean cmean./smean tdis_prb(cmean./smean,nobs) lbounds ubounds  ];
end;


fprintf(fid,'\n');
fprintf(fid,'Spatial autoregressive Model Estimates \n');
if (nflag == 1)
fprintf(fid,'Dependent Variable = %16s \n',vnames(1,:));
end;
fprintf(fid,'R-squared          = %9.4f \n',results.rsqr);
fprintf(fid,'Rbar-squared       = %9.4f \n',results.rbar);
fprintf(fid,'sigma^2            = %9.4f \n',results.sig2_v);
fprintf(fid,'Nobs, Nvars        = %6d,%6d \n',results.nobs,results.nvar);
fprintf(fid,'log-likelihood     = %16.8g \n',results.lik);
%fprintf(fid,'# of iterations    = %6d   \n',results.iter);
fprintf(fid,'min and max rho    = %9.4f,%9.4f \n',results.rmin,results.rmax);
% print timing information
fprintf(fid,'total time in secs = %9.4f \n',results.time);
if results.time1 ~= 0
fprintf(fid,'time for lndet     = %9.4f \n',results.time1);
end;
if results.time2 ~= 0
fprintf(fid,'time for eigs      = %9.4f \n',results.time2);
end;
if results.time3 ~= 0
fprintf(fid,'time for t-stats   = %9.4f \n',results.time3);
end;
if results.time5 ~= 0
fprintf(fid,'time for x-impacts = %9.4f \n',results.time5);
fprintf(fid,'# draws  x-impacts = %9d   \n',results.ndraw);
end;    

if results.lflag == 0
fprintf(fid,'No lndet approximation used \n');
end;
% put in information regarding Pace and Barry approximations
if results.lflag == 1
fprintf(fid,'Pace and Barry, 1999 MC lndet approximation used \n');
fprintf(fid,'order for MC appr  = %6d  \n',results.order);
fprintf(fid,'iter  for MC appr  = %6d  \n',results.miter);
end;
if results.lflag == 2
fprintf(fid,'Pace and Barry, 1998 spline lndet approximation used \n');
end;

fprintf(fid,'***************************************************************\n');

bout = [results.beta
        results.rho];
    
% now print coefficient estimates, t-statistics and probabilities
tout = norm_prb(results.tstat); % find asymptotic z (normal) probabilities
tmp = [bout results.tstat tout];  % matrix to be printed
% column labels for printing results
bstring = 'Coefficient'; tstring = 'Asymptot t-stat'; pstring = 'z-probability';
cnames = strvcat(bstring,tstring,pstring);
in.cnames = cnames;
in.rnames = Vname;
in.fmt = '%16.6f';
in.fid = fid;
mprint(tmp,in);

% now print x-effects estimates

bstring = 'Coefficient'; 
tstring = 't-stat'; 
pstring = 't-prob';
lstring = 'lower 01';
ustring = 'upper 99';
cnames = strvcat(bstring,tstring,pstring,lstring,ustring);
ini.cnames = cnames;
ini.width = 2000;

% print effects estimates
if cflag == 1
vnameso = strvcat(Vname(3:end-1,:));
elseif cflag == 0
vnameso = strvcat(Vname(2:end-1,:));    
end
ini.rnames = strvcat('Direct  ',vnameso);
ini.fmt = '%16.6f';
ini.fid = fid;

% set up print out matrix
printout = direct_out;
mprint(printout,ini);

printout = indirect_out;
ini.rnames = strvcat('Indirect',vnameso);
mprint(printout,ini);

printout = total_out;
ini.rnames = strvcat('Total   ',vnameso);
mprint(printout,ini);




        
% <=================== end of sar case




otherwise
error('results structure not known by prt_sar function');
end;





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

% written by:
% James P. LeSage, 3/2010
% Dept of Finance & Economics
% Texas State University-San Marcos
% 601 University Drive
% San Marcos, TX 78666
% jlesage@spatial-econometrics.com


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


























