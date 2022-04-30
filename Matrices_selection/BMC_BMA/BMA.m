% PURPOSE: An example of model comparison using sar_c() function to compare various weight matrix specifications
%          (on a small data set)                  
%---------------------------------------------------
% USAGE: compare_weights2
%---------------------------------------------------


clear all;

load elect.dat;                    % load data on votes
latt = elect(:,5);
long = elect(:,6);
n = length(latt);

% create W-matrix based on nearest 3 neighbors
W3 = make_neighborsw(latt,long,3);


% generate an sar model based on 4 nearest neighbors
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


% estimate 5 models using W1 to W5 as weight matrices

% run 5 homoscedastic models
ndraw = 2500;
nomit = 500;
prior.thin=1;
W1 = make_neighborsw(latt,long,1); % create W-matrix based on nearest 1 neighbor
results1 = sar_g(y,x,W1,ndraw,nomit,prior);
prt(results1,vnames)
% prt(results1,vnames);  NOTE: we cannot use prt here since there are no estimates to print           
W2 = make_neighborsw(latt,long,2); % create W-matrix based on nearest 2 neighbors
results2 = sar_g(y,x,W2,ndraw,nomit,prior);
prt(results2,vnames)
% prt(results2,vnames);
W3 = make_neighborsw(latt,long,3); % create W-matrix based on nearest 3 neighbors
results3 = sar_g(y,x,W3,ndraw,nomit,prior);
prt(results3,vnames)
% prt(results3,vnames);
W4 = make_neighborsw(latt,long,4); % create W-matrix based on nearest 4 neighbors
results4 = sar_g(y,x,W4,ndraw,nomit,prior);
prt(results4,vnames)
% prt(results4,vnames);
W5 = make_neighborsw(latt,long,5); % create W-matrix based on nearest 5 neighbors
results5 = sar_g(y,x,W5,ndraw,nomit,prior);
prt(results5,vnames)

W=[W1 W2 W3 W4 W5];
res=sar_bma_g(y,x,W,ndraw,nomit);

% prt(results5,vnames);







