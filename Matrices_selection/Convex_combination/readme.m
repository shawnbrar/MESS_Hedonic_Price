% Readme file

% This directory contains all Matlab codes needed to estimate a SAR model with a convex 
% combination of connectivity matrices as well as the BMA estimates. 

demo.m                  % Demo file with all functions
model_probs.m :         % function used to compute posterior model probabilities used in the Bayesian Model Averaging procedure.
prt_sar_conv_bma_g.m :  % Printing function of the BMA estimates.
prt_sar_conv_g.m:       % Printing function for the SAR model with a convex combination of W matrices estimation procedure.
prt_sar_g.m:            % Printing function for the conventional SAR model estimated by MCMC, modified so that log-marginal
                        % likelihood is constructed in the same way as in the SAR model with convex combination of W matrices.
sar_conv_bma_g.m :      % Estimation  of the Bayesian Moving Averaged estimates.
sar_conv_g.m :          % Estimation of the SAR model with a convex combination of W matrices
sar_g.m :               % Estimation of the conventional SAR model by MCMC, modified so that log-marginal
                        % likelihood is constructed in the same way as in the SAR model with convex combination of W matrices.