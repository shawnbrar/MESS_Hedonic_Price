% file: cross_section_demo3.m, demonstrates comparison of weight matrices
clear all;
rng(86573);
load uscounties.data;
% a matrix now exist named uscounties
% the matrix contains 11 columns of county-level data
% col 1  FIPS     
% col 2  LATITUDE 
% col 3  LONGITUDE
% col 4  POP1990  
% col 5  1987_PCI (per capita income)
% col 6  1988_PCI 
% col 7  1989_PCI 
% col 8  1990_PCI 
% col 9 1991_PCI 
% col 10 1992_PCI 
% col 11 1993_PCI 
[n,k] = size(uscounties); % find the size of the matrix
pci1987 = uscounties(:,5);  % extract the 5th column from the data matrix 
pci1993 = uscounties(:,11); % creates an n x 1 column vector 
pop1990 = uscounties(:,4);
% calculate growth rate of per capita income over the 1987 to 1993 period
pci_growth = log(pci1993) - log(pci1987);
% make these annualized growth rates
pci_growth = pci_growth/7;
% do a growth regression
% which involves regressing the growth rate on the (logged) initial level
xmatrix = [log(pci1987) log(pop1990)];
% run SDM model
latt = uscounties(:,2); % extract latt-long coordinates
long = uscounties(:,3);
W(1).matrix = zeros(n,n);
neigh = [];
for ii=4:16
    W(ii-3).matrix = make_neighborsw(latt,long,ii);
    neigh = [neigh
             ii];
end

lmarginal = [];
nweights = size(neigh,1);

for i=1:nweights

res(i).result = lmarginal_cross_section(pci_growth,xmatrix,W(i).matrix);

lmarginal = [lmarginal
             res(i).result.lmarginal];

end

probs = model_probs(lmarginal);

probs_matrix = reshape(probs,nweights,3);

in.fmt = '%16.4f';
in.cnames = strvcat('slx','sdm','sdem');
in.rnames = strvcat('#neighbors',num2str(neigh));
mprint(probs_matrix,in);

% run the best model
% 1st we need to figure out which model is best
% and how many neighbors to use
[nprob,nmax] = max(probs_matrix,[],1);

[prob,mmax] = max(nprob);

if mmax == 1
    model = 'slx';
elseif mmax == 2
    model = 'sdm';
elseif mmax == 3
    model = 'sdem';
end

neighbors_index = nmax(mmax);
    
num_neighbors = neigh(nmax(mmax));


W = make_neighborsw(latt,long,num_neighbors);
ndraw = 2500;
nomit = 500;
prior.novi_flag = 1;
prior.model = 0;
T = 1;

switch model
    
    case{'slx'}        

result = slx_panel_FE_g(pci_growth,[ones(n,1) xmatrix],W,T,ndraw,nomit,prior);

    case{'sdm'}        

result = sdm_panel_FE_g(pci_growth,[ones(n,1) xmatrix],W,T,ndraw,nomit,prior);

    case{'sdem'}        

result = sdem_panel_FE_g(pci_growth,[ones(n,1) xmatrix],W,T,ndraw,nomit,prior);

    otherwise
        % do nothing
       disp('Unknown model'); 

end

vnames = strvcat('income growth','constant','log(pci0)','log(pop0)');
prt_panel(result,vnames);

