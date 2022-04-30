%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAR model in cross section%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
load data.mat
n = size(x)(1);
cons=ones(n, 1);
x=[cons x];
x = x(:, [1:12, 14]); % deleting the LOTZ variable because of multicollinearity
load spt_mat;
W = spt_mat;
Wi = normw(W); % row normalize W
clear spt_mat;

info.lflag=0;
vnames=strvcat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MESS Estimation %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ws = W/eigs(W)(1, 1);
cd '~/Documents/DASEE_Study/Sem 3/Adv_Tp_Spatial_Stats/Octave_Code/MESS'
res1b=mess_nls(y, x, Ws); % beta, alpha

% Impacts %%%%%%%%%%%%%%%%%%%%
setappdata(0, 'Ws', Ws);
setappdata(0, 'alpha', res1b.alpha);
use_betas = num2cell(res1b.beta(2:end))

function k = Impactz(betaz)
  globalWs = getappdata(0, 'Ws');
  globalalpha = getappdata(0, 'alpha');
  n = size(globalWs)(1, 1);
  SrW = expm(-globalalpha * globalWs) * betaz;
  Ad = (1/n)*trace(SrW);       % direct impacts
  At = (1/n)*sum(sum(SrW, 2)); % total impacts
  Ai = At - Ad;                % indirect impacts
  k = [Ad Ai At];
endfunction

impacts = cell2mat(cellfun(@Impactz, use_betas, 'UniformOutput', false))

%% Inference on Impacts

cov = res1b.cov;
ex = expm(-res1b.alpha*W);
Deriv = zeros(size(x, 2) + 1, 1);
Deriv(2, 1) = (1/n)*trace(ex);
Deriv(end, 1) = (1/n)*trace((-Ws)*ex)*res1b.beta(2);
den = (transpose(Deriv) * cov * Deriv)^0.5;
tstat = impacts(1, 1)/den;

cov = res1b.cov;
ex = expm(-res1b.alpha*Ws);
covari = size(res1b.beta(2:end))(1);
%%%% Std dev of Direct Impacts
deriv = [eye(covari);res1b.beta(2:end)'];
trace1 = 1/n*trace(ex);
trace2 = 1/n*trace(-Ws*ex);
deriv = [repmat(trace1 , covari, 1);trace2].*deriv;
denom_vector = zeros(12, 1);
for i = 1:covari
  % std dev for the direct impacts
  denom_vector(i, 1) = (transpose(deriv(:, i)) * cov(2:end, 2:end) * deriv(:, i))^(0.5);
endfor

%%%% Std dev of Total Impacts
deriv2 = [eye(covari);res1b.beta(2:end)'];
ex2 = ex - diag(diag(ex));
trace1 = 1/n*sum(sum(ex));
trace2 = 1/n*sum(sum(-Ws*ex));
deriv = [repmat(trace1 , covari, 1);trace2].*deriv;
denom_vector2 = zeros(12, 1)
for i = 1:covari
  % std dev for the total impacts
  denom_vector2(i, 1) = (transpose(deriv2(:, i)) * cov(2:end, 2:end) * deriv2(:, i))^(0.5);
endfor

%%%% Std dev of Indirect Impacts
deriv3 = [eye(covari);res1b.beta(2:end)'];
ex2 = ex - diag(diag(ex));
trace1 = 1/n*sum(sum(ex2));
trace2 = 1/n*sum(sum(-Ws*ex2));
deriv3 = [repmat(trace1 , covari, 1);trace2].*deriv3;
denom_vector3 = zeros(12, 1)
for i = 1:covari
  % std dev for the indirect impacts
  denom_vector3(i, 1) = (transpose(deriv3(:, i)) * cov(2:end, 2:end) * deriv3(:, i))^(0.5);
endfor

std_dev = [denom_vector denom_vector3 denom_vector2]; % std dev value of direct, indirect and total effects 
tvals = impacts./std_dev; % t value of direct, indirect and total effects 
pvals = norm_prb(tvals); % significane probability of direct, indirect and total effects
clear deriv deriv2 deriv3 trace1 trace2 ex ex2 denom_vector denom_vector2 denom_vector3;
pvals; % p value for significance of direct, indirect and total effects 

model_coef = [res1b.beta;res1b.alpha]
model_mat = [model_coef res1b.tstat norm_prb(res1b.tstat)]