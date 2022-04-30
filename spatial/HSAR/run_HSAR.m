%Paper of Aquaro, Bailey and Pearan (JAE 2020): Estimation and inference for spatial models with heterogeneous coefficients: An application to US house prices
%Data and program come from Michele Aquaro

% Estimate HSAR models for quarterly real house price changes in the USA 
% at Metropolitan Statistical Areas (MSAs) over the period 1975:Q1?2014:Q4. 
% Use of 338 MSA over 159 quarters

% Account for the presence of common factor that might affect  house price
% changes

%\pi_{it}= a_i + \Psi_{0,i}\sum_{j=1}^N w_{ij}\pi_{jt} + \Psi_{1,i}\sum_{j=1}^N w_{ij}\pi_{j,t-1} + \lambda_i \pi_{i,t-1}
%+ \beta_i^{pop}gpop_{it} + \beta_i^{inc}ginc_{i,t} + \varepsilon_{it}

% Changes in house prices depend on current house prices changes in the
% neighborhood, past price changes in the neighborhood, past price changes
% in the MSA, on percent rate of change of population and rate of change of (real)income per capita

%% Choice of the W matrix
%Geodesic distance between each pair of latitude/longitude coordinates for the MSAs included in our sample 
%These coordinates correspond to the center of the polygon implied by each MSA. 
% Then, we determine a specific radius threshold, d (miles), within which MSAs are considered to be neighbors. In this case, the relevant entries in the un-normalized weight matrix W0 are set to unity. The MSAs that fall outside this radius are labeled non-neighbors and their corresponding entries in W0 are set to zero. Finally, we row-normalize W0 and obtain W, which is used in Equation 30.
%3 versions of W constructed with the radius threshold values of d = 75, 100 and 125 miles. 

% Here, estimation with 75 miles

load data_75.mat
%a_x = ones(N, T - 1, K); %intercept
%a_x(:, :, 2) = m_pp(:, 2:T);% population change (in percent)
%a_x(:, :, 3) = m_ic(:, 2:T); % real income per capita change (in percent)
%a_x(:, :, 4) = m_hs(:, 1:(T - 1)); % Wy1 (spatio-temporal y)
%a_x(:, :, 5) = m_hp(:, 1:(T - 1)); % y1 (lag of y)

% Estimate the HSAR model
%% upper bounds
[N,T,K]=size(a_x);
m_ub_x = NaN(N, K); 
m_ub_x(:, 1) = Inf(N, 1); % intercept
m_ub_x(:, 2) = Inf(N, 1); % pp
m_ub_x(:, 3) = Inf(N, 1); % ic
m_ub_x(:, 4) = 0.995 * ones(N, 1); % Wy1
m_ub_x(:, 5) = 0.995 * ones(N, 1); % y1

%% lower bounds
m_lb_x = -m_ub_x;

% define parameter space for the psi's and sgmsq's
v_ub_Wy    = 0.995 * ones(N, 1);
v_lb_Wy    = -v_ub_Wy;
v_ub_sgmsq = Inf(N, 1);
v_lb_sgmsq = zeros(N, 1);

% put all bounds together into a single array of order (N,K+2,2)
m_lb = [v_lb_Wy m_lb_x v_lb_sgmsq];
m_ub = [v_ub_Wy m_ub_x v_ub_sgmsq];
a_b = NaN(N, K + 2, 2);
a_b(:, : , 1) = m_lb;
a_b(:, : , 2) = m_ub;

% -------------------------------------------------------------------------- %
% initial values for the optimisation procedure
v_theta_ini_tr = zeros(1, K + 2); % (1,K+2)
v_theta_ini_tr(end) = 1; % sgmsq
m_theta_ini = repmat(v_theta_ini_tr, [N 1]); % (N,K+2)

% -------------------------------------------------------------------------- %
% estimate HSAR
results = hsar(m_hp, a_x, m_W, a_b, m_theta_ini);
% generate residuals
[N T] = size(m_hp);
m_e = zeros(N, T);
for ii = 1:N
     psi_hat = results.m_theta(ii, 1); % (1,1)
     v_beta_hat = results.m_theta(ii, 2:(end - 1))'; % (K,1)
     for t = 1:T
          y = m_hp(ii, t);
          v_x = squeeze(a_x(ii, t, :)); % (K,1)
          y_hat = v_beta_hat' * v_x;
          m_e(ii, t) = y - y_hat;
     end
end

% -------------------------------------------------------------------------- %
% save results
save(sprintf('estimates_W%03d.mat', miles), ...
     'v_id', ...
     'v_reg', ...
     'results', ...
     'v_time', ... included for completeness, not needed
     'm_e', ...
     'm_W' ...
     )
%  
%  results.theta(:,1) psi
%  results.theta(:,2) alpha_i
%  results.theta(:,3) pop
%  results.theta(:,4)=inc
%  results.theta(:,5)= Wy1
%  results.theta(:,6)=y1;
 
 


