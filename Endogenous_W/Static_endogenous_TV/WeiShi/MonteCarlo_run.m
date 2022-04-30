% 2015/10/22

clearvars
clc
global n T theta_select Ry0 R SNR

R               = 10;     % number of Monte Carlo iterations
theta_select    = 1;      % parameter choice, 1,2,3,4
Ry0             = 2;      % true number of factors in the y equation
SNR             = 1;      % variance of v and e, 1,2,4,9

n               = 49;
T               = 25;
run('MonteCarlo.m')