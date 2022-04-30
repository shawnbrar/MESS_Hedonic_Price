% PURPOSE: An example of model comparison using sar_c() function to compare various weight matrix specifications
%          on a homoscedastic sar model   
% (see compare_weights2 for non-homoscedastic models)               
%---------------------------------------------------
% USAGE: compare_weights
%---------------------------------------------------

clear all;

load elect.dat;                    % load data on votes
latt = elect(:,5);
long = elect(:,6);
n = length(latt);
% x = [ones(n,1) x1 x2 x3];
% clear x1; clear x2; clear x3;
clear elect;                % conserve on RAM memory


% create W-matrix based on nearest 3 neighbors
W3 = make_neighborsw(latt,long,3);

% generate an sar model based on 4 nearest neighbors
n = length(latt);
IN = eye(n); 
rho = 0.7;  % true value of rho
sige = 0.2;
k = 3;
x = randn(n,k);
x(:,1) = ones(n,1);
beta(1,1) = -1.0;
beta(2,1) = 0.0;
beta(3,1) = 1.0;

vnames = strvcat('y','constant','x1','x2');
    
% sar model generated here
% based on nearest 3-neighbors W-matrix, (W3 from above)

y = (IN-rho*W3)\(x*beta) + (IN-rho*W3)\(randn(n,1)*sqrt(sige)); 


% compute log-marginal posteriors for 5
% homoscedastic models using W1 to W5 as weight matrices

% run 5 homoscedastic models
ndraw = 2500;
nomit = 500;
prior.lflag=0;

W1 = make_neighborsw(latt,long,1); % create W-matrix based on nearest 1 neighbor
results1 = sar_g(y,x,W1,ndraw,nomit,prior);
% prt(results1,vnames);  NOTE: we cannot use prt here since there are no estimates to print           
W2 = make_neighborsw(latt,long,2); % create W-matrix based on nearest 2 neighbors
results2 = sar_g(y,x,W2,ndraw,nomit,prior);
% prt(results2,vnames);
results3 = sar_g(y,x,W3,ndraw,nomit,prior);
% prt(results3,vnames);
W4 = make_neighborsw(latt,long,4); % create W-matrix based on nearest 4 neighbors
results4 = sar_g(y,x,W4,ndraw,nomit,prior);
% prt(results4,vnames);
W5 = make_neighborsw(latt,long,5); % create W-matrix based on nearest 5 neighbors
results5 = sar_g(y,x,W5,ndraw,nomit,prior);
% prt(results5,vnames);

% compare 5 homoscedastic models based on 5 different weight matrices
fprintf(1,'posterior probabilities for 5 homooscedastic models \n');
fprintf(1,'based on W-matrices for neighbors 1 to 5 \n');
lmarginal = [results1.logm_sar results2.logm_sar results3.logm_sar results4.logm_sar results5.logm_sar];

adj = max(lmarginal);
madj = lmarginal - adj;

xx = exp(madj);

% compute posterior probabilities
psum = sum(xx);
probs = xx/psum;
 probs=probs';
 lmar=lmarginal';
%  
% results.probs = probs';
% results.lmarginal = lmarginal';
% probs = model_probs(results1, results2, results3, results4, results5);
out=[probs,lmar];
rnames = strvcat('Homoscedastic Models');
for j=1:5
    rnames = strvcat(rnames,['neighbors ',num2str(j)]);
end;
in.rnames = rnames;
in.cnames = strvcat('Model Probabilities','log-marginal');
fprintf(1,'True model is based on 3 neighbors \n');
mprint(out,in);

