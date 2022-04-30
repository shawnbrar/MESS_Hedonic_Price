% model_comparison_chapter5p2.m demo file
clear all;
[cclaims,b] = xlsread('../demo_data/weekly.xlsx',1);
% read data from sheet 1 of Excel spreadsheet
% growth rate of unemployment 2019-2020 from same week, previous year
snames = strvcat(b(2:end,1)); % 48 state names
tnames = strvcat(b(1,2:end)); % 51 week labels
[N,T] = size(cclaims);
[jobposts,b] = xlsread('../demo_data/weekly.xlsx',2);
% read data from sheet 2 of Excel spreadsheet
% change in job offers from 1st week of 2020
[athome,b] = xlsread('../demo_data/weekly.xlsx',3);
% read data from sheet 3 of Excel spreadsheet
% growth rate of percent population at home
% 2019-2020 from same week, previous year
[a,b] = xlsread('../demo_data/Wcont48.xlsx');
% 48 x 48 contiguity matrix for states
Wc = normw(a);
% miles of borders in common
[a,b] = xlsread('../demo_data/states_borders.xlsx');
Wmiles = a(:,2:end);
% only upper triangular
% so we make it symmetric
for i=1:48
    for j=1:48
        if Wmiles(i,j) > 0
            Wmiles(j,i) = Wmiles(i,j);
        end
    end
end

Wb = normw(Wmiles);

y = vec(cclaims);
x = [vec(jobposts) vec(athome)];

model = 3; % both state and time fixed effects
[ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,x,N,T,model);

info.lflag = 0; % exact log-determinant
% info.lflag = 1; uses Pace and Barry approximation
% which is faster for large data samples

result1c = lmarginal_static_panel(ywith,xwith,Wc,N,T,info); 
result1b = lmarginal_static_panel(ywith,xwith,Wb,N,T,info); 

fprintf(1,'state and time fixed effects: log marginal likelihoods and model probabilities \n');
in.cnames = strvcat('Wc log-marginal','model probs');
in.rnames = strvcat('model','slx','sdm','sdem');
in.width = 10000;
in.fmt = '%10.4f';
out1 = [result1c.lmarginal result1c.probs];
mprint(out1,in);

fprintf(1,'state and time fixed effects: log marginal likelihoods and model probabilities \n');
in.cnames = strvcat('Wb log-marginal','model probs');
in.rnames = strvcat('model','slx','sdm','sdem');
in.width = 10000;
in.fmt = '%10.4f';
out2 = [result1b.lmarginal result1b.probs];
mprint(out2,in);

% compare all models
lmarginals = [result1c.lmarginal
              result1b.lmarginal];
          
probs = model_probs(lmarginals);
% rearrange for pretty printing
mprobs = reshape(probs,3,2);

out = [result1c.lmarginal mprobs(:,1) result1b.lmarginal mprobs(:,2)];
    
fprintf(1,'Comparison of Wc and Wb \n');
in.cnames = strvcat('Wc log-marginal','Wc model probs','Wb log-marginal','Wb model probs');
in.rnames = strvcat('model','slx','sdm','sdem');
in.width = 10000;
in.fmt = '%10.4f';
mprint(out,in);

% % =============================================
% % estimate SDM models based on Wc and Wb
% 
ndraw = 3000;
nomit = 1000;

prior.model = 3;
prior.novi_flag = 1;
result1 = sdm_panel_FE_g(y,x,Wc,T,ndraw,nomit,prior);
vnames = strvcat('cclaims','jobposts','athome');
prt_panel(result1,vnames);
% 
result2 = sdm_panel_FE_g(y,x,Wb,T,ndraw,nomit,prior);
vnames = strvcat('cclaims','jobposts','athome');
prt_panel(result2,vnames);

prior2.model = 3;
prior2.rval = 4;
result1 = sdm_panel_FE_g(y,x,Wb,T,ndraw,nomit,prior2);
vnames = strvcat('cclaims','jobposts','athome');
prt_panel(result1,vnames);
% 
result2 = sdem_panel_FE_g(y,x,Wb,T,ndraw,nomit,prior2);
vnames = strvcat('cclaims','jobposts','athome');
prt_panel(result2,vnames);
