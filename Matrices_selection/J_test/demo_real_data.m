clear all
load data59.txt;
data59(13,:)=[];
n = length(data59);
xlong = data59(:,1);
ylat = data59(:,2);
lny = data59(:,3);
lns = data59(:,4);
lnngd = data59(:,5);
lnn = log(exp(lnngd)-0.05);
lnrd = data59(:,6);
%Old human capital
H = data59(:,7);
h = data59(:,8);
%new human capital measure, based on Cohen data based in PSE
load Human_cap_Cohen_database_PSE %Without Korea
H_90=H_1990; %Human capital in 1990
H_00=H_2000; %Human capital in 2000

info.lflag = 0;
%% Creation of interaction matrices
% Genetic Distance
load DistGen3.txt;

i=DistGen3(:,1);
j=DistGen3(:,2);
fstw=DistGen3(:,6);
% building of W
A=[i j fstw];
S=spconvert(A);
B=full(S);
dgW=B+B';


% Final specification
dgWinv=exp(-1*(dgW/100))-eye(59);


% Remove the 13�me observation from the sample because in the linguistic interaction scheme,
% it does not interact with anyone (line full of zeros).
dgWinv(13,:)=[];
dgWinv(:,13)=[];
% gen1=diag(h)*normw(dgWinv);
gen1=diag(H_90)*normw(dgWinv);
% gen1=diag(H_00)*normw(dgWinv);


C3=max(sum(gen1));
R3=max(sum(gen1'));
gen2=gen1./min(R3,C3);

clear gen1 R2 C3A S B dgW dgWinv fstw i j 
%%Linguistic proximity
load linguis.txt;
% Remove the 13�me observation from the sample because it does not interact
% with anyone (line full of zeros).

linguis(13,:)=[];
linguis(:,13)=[];
% ling = diag(h)*normw(linguis);
ling = diag(H_90)*normw(linguis);
% ling = diag(H_00)*normw(linguis);


Cl=max(sum(ling));
Rl=max(sum(ling'));
lingui=ling./min(Cl,Rl);
clear ling linguis Cl Rl 


%% Trade matrix (unconstrained) 
load trade.txt; % Unconstrained trade matrix 
% Remove the 13�me observation from the sample because in the linguistic interaction scheme,
% it does not interact with anyone (line full of zeros).
trade(13,:)=[];
trade(:,13)=[];
% tradee = diag(h)*normw(trade);
tradee = diag(H_90)*normw(trade);
%  tradee = diag(H_00)*normw(trade);


Ct=max(sum(tradee));
Rt=max(sum(tradee'));
tradeT=tradee./min(Ct,Rt);
clear tradee Ct Rt


x = [ones(n,1) lnrd lnn 0.5*(lns-lnngd)];
vnames = strvcat('lnTFP95','constant','lnrd','lnn','lnsngd');
% Estimation for the interaction matrix based on geneoalogical similarity.
Wgx=gen2*x;
Wg2x=gen2*Wgx;
Q=[Wgx Wg2x x];
P1=gen2;

% tfpgmmegen2=gmmestimation_sar_he(lny-0.5*(lns-lnngd),xsolowAHTFP,gen2,Q,P1);
tfpgmmegen2=gmmestimation_sar_he_old(lny-0.5*(lns-lnngd),x,gen2,Q,P1);


% Estimation for an interaction matrix based on trade flows
% solowAHTFPtradeT = sarq(lny-0.5*(lns-lnngd),xsolowAHTFP,tradeT,info);
% prt(solowAHTFPtradeT,vsolowAHTFP);
Wttx=tradeT*x;
Wtt2x=tradeT*Wttx;
Q=[Wttx Wtt2x x];
P1=tradeT;
% tfpgmmetradeT=gmmestimation_sar_he(lny-0.5*(lns-lnngd),xsolowAHTFP,tradeT,Q,P1);
tfpgmmetradeT=gmmestimation_sar_he_old(lny-0.5*(lns-lnngd),x,tradeT,Q,P1);

% Estimation with interaction matrix modeled by linguistic proximity
   Wlx=lingui*x;
    Wl2x=lingui*Wlx;
    Q=[Wlx Wl2x x];
    P1=lingui;
%     tfpgmmelingui=gmmestimation_sar_he(lny-0.5*(lns-lnngd),xsolowAHTFP,lingui,Q,P1);
    tfpgmmelingui=gmmestimation_sar_he_old(lny-0.5*(lns-lnngd),x,lingui,Q,P1);
    
    display('Hageman and J-test based on the unconstrained TFP model');
   
   W1=tradeT;
W2=gen2; % Genetic
W3=lingui; % linguistic
res_1=tfpgmmetradeT;
res_2=tfpgmmegen2;
res_3=tfpgmmelingui;

B=99;
n=58;
y=lny-0.5*(lns-lnngd);
[stat1_p1, stat2_p1,stat3_p1,signi_p1]=Jp1_hetero_3(y,W1,W2,W3,x,res_1,res_2,res_3,B);

[val_p1, hag_p1]=min([stat1_p1 stat2_p1 stat3_p1],[],2);
% [val_p2, hag_p2]=min([stat1_p2 stat2_p2 stat3_p2],[],2);
%  [val_p3, hag_p3]=min([stat1_p3 stat2_p3 stat3_p3],[],2);
choice=[1 'trade     ';
        2 'genetic   ';
        3,'linguistic']
    
display('Results concerning J test based on RGMM')
display('for first predictor')

out11=fprintf('Hagemann procedure selects matrix %s.\n',choice(hag_p1,2:end));
out12=fprintf('Associated value of J test is %.3f.\n',val_p1);
out12a=fprintf('Critical value for J-test of the selected matrix using bootstraped distribution is %2.3f.\n',signi_p1.qu_t);
out12b1=fprintf(' value of J test for first alternative model is %.3f.\n',stat2_p1);
out12b2=fprintf('Critical value for J-test of the first alternative model using bootstraped distribution is %2.3f.\n',signi_p1.pval_M1alt);
out12c1=fprintf(' value of J test for second alternative model is %.3f.\n',stat3_p1);
out12c2=fprintf('Critical value for J-test of the second alternative model using bootstraped distribution is %2.3f.\n',signi_p1.pval_M2alt);
out13=fprintf('Pseudo p_value associated to the J test using %d replications is %.3f.\n',B,signi_p1.pval_J);
out14=fprintf('P-value based on asymptotic distribution (chi squqare with 2 dof) is %.3f.\n',1-chis_prb(val_p1,2));
out15=fprintf('J-test (predictor 1) for W1 (W1 true) and associated p-value (asymptotic) is %.3f    %.3f.\n',stat1_p1,1-chis_prb(stat1_p1,2));
% out15a=fprintf('J-test (predictor 1) for W1 (W1 true) and associated bootstraped p-value is %.3f with a p-value of   %.3f.\n',stat1_p1_b,signi_p1_b.pva1);
out16=fprintf('J-test (predictor 1) for W2 (W2 true) and associated p-value (asymptotic) is %.3f    %.3f.\n',stat2_p1,1-chis_prb(stat2_p1,2));
% out16a=fprintf('J-test (predictor 1) for W2 (W1 true) and associated bootstraped p-value is %.3f with a p-value of   %.3f.\n',stat2_p1_b,signi_p1_b.pva2);
out17=fprintf('J-test (predictor 1) for W3 (W3 true) and associated p-value (asymptotic) is %.3f    %.3f.\n',stat3_p1,1-chis_prb(stat3_p1,2));

%[stat12, stat22,stat32,signi2]=Jp2_hetero_3(y,W1,W2,W3,x,res_1,res_2,res_3,B)

