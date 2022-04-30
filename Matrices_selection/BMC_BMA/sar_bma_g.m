function results=sar_bma_g(y,x,W,ndraw,nomit,prior)
  % Computes Bayesian Model Averaged estimated of the SAR specification. 
  % W includes all matrices, concatenated by rows : W=[W1 W2 W3];

if ndraw <= 1000
       error('sar_bma_g: ndraw<=1000, increase ndraw to at least 10000');
end

  
  [n,k]=size(x);
if nargin == 5
    results.thin=1;
    end
  
results.meth="sar_bma_g";
results.nobs  = n;
results.nvar  = k;
results.y = y; 
% check if the user handled the intercept term okay
    if sum(x(:,1)) ~= n
        tst = sum(x); % we may have no intercept term
        ind = find(tst == n); % we do have an intercept term
        if length(ind) > 0
            error('sar_bma_g: intercept term must be in first column of the x-matrix');
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
    
[n1,m]=size(W);
  np=m/n1;
  if np~= round(np)
        error('sar_bma_g: wrong sized W-matrices');
    elseif n1 ~= n
        error('sar_bma_g: wrong sized W-matrices');
    elseif np < 2
        error('sar_bma_g: only one W-matrix');
    end
   results.np=np; 
  for i =1:np
    begi=(i-1)*n+1;
    endi=i*n;
    Wmatrix{i}=W(:,begi:endi);
    end
   for iter = 1:np
    res = sar_g(y,x,Wmatrix{iter},ndraw,nomit);
    
    temp(iter) = res;
 end
 results.rmin=temp(1).rmin;
  results.rmax=temp(1).rmax;
 
 results.thin=temp(1).thin;;
 thin=results.thin;
 if thin ~= 1
    neff = (ndraw-nomit)/prior.thin;
 else
     neff = ndraw-nomit;
 end
 bsave = zeros(np,neff,k);

psave = zeros(np,neff);
ssave = zeros(np,neff);
smean=zeros(np,1);
dsave = zeros(np,neff,pp);
isave = zeros(np,neff,pp);
tsave = zeros(np,neff,pp);
rho_mean = zeros(np,1);
beta_mean = zeros(np,k);
logmarginal = zeros(np,1);


for ii=1:np
    
    
    logmarginal(ii,1) = temp(ii).logmarginal;
    rho_mean(ii,1) = temp(ii).rho;
    beta_mean(ii,:)=temp(ii).beta';
    dsave(ii,:,:) = temp(ii).direct;
    isave(ii,:,:) = temp(ii).indirect;
    tsave(ii,:,:)=temp(ii).total;
    bsave(ii,:,:) = temp(ii).bdraw;
    psave(ii,:) = temp(ii).pdraw';
    ssave(ii,:) = temp(ii).sdraw';
    smean(ii,1)=temp(ii).sigma;
    ttime(ii,1)=temp(ii).time;
       
end 

% Computation of posterior probabilities
adj = max(max(logmarginal));
madj = logmarginal - adj;

xx = exp(madj);

% compute posterior probabilities
psum = sum(xx);
probs = xx/psum;

prt_out3 = [logmarginal probs rho_mean];

cnames = strvcat('logmarginal','Prob','rho');


in.cnames = cnames;

% rnames = strvcat('Models');
rnames = strvcat('Models');
for i=1:np
    rnames = strvcat(rnames,['W' num2str(i)]); 
end


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

bma_logm = sum(probs.*logmarginal);
prt_out3 = [prt_out3
            bma_logm  1 bma_rho 
            logmarginal(hi,1) probs(hi,:) rho_mean(hi,:)];
 in.rnames = strvcat(rnames, 'BMA','highest');
mprint(prt_out3,in);


%% Computation of the BMA values for coefficients and effects 

    bma_beta=(probs'*beta_mean)';
    
    
results.probs=probs;
results.rho=bma_rho;
results.beta=bma_beta;
results.sige = bma_sige;


bma_b=zeros(neff,k);
bma_dir=zeros(neff,pp);
bma_ind=zeros(neff,pp);
bma_tot=zeros(neff,pp);
bma_r=zeros(neff,1);
bma_s=zeros(neff,1);


for j=1:k
    bma_b(:,j)=sum(matmul(probs,bsave(:,:,j)));
end
    bma_r(:,1)=sum(matmul(probs,psave));
    bma_s(:,1)=sum(matmul(probs,ssave));

for j = 1:pp
    bma_dir(:,j)=sum(matmul(probs,dsave(:,:,j)));
    bma_ind(:,j)=sum(matmul(probs,isave(:,:,j)));
    bma_tot(:,j)=sum(matmul(probs,tsave(:,:,j)));
end


results.total=bma_tot;
results.direct=bma_dir;
results.indirect=bma_ind;
results.meth= 'sar_bma_g';
results.tflag = 'plevel';
results.ndraw=ndraw;
results.nomit=nomit;

results.bdraw=bma_b;
results.pdraw=bma_r;
results.sdraw=bma_s;
results.lmarginal=bma_logm;

end