% Kalman Filter
% Input : A, B, C - System Matrices
% Input : V, W - Process and Measurement Noise Covariance
% Input : x_bar_est : state estimate, y - observation vector
% Input : u - previous input
% Output : x_upd - State Estimate
% Output : P_upd - Posterior State Covariance Matrix
function [x_upd, P_upd] = KF_prd(A,B,C,V,W,P_prd,x_prd,y,u)
    n = length(x_prd);
    % x_prd = A*x_est; % State Prediction
    % P_prd = A*P*(A.') + V; % Prior Covariance
    K = P_prd*(C.')/(C*P_prd*(C.')+W); % Kalman Gain
    x_est = x_prd + K*(y-C*x_prd); % State Estimate
    P_est = (eye(n)-K*C)*P_prd; % Posterior Covariance
    x_upd = A*x_est + B*u; % Next state prediction
    P_upd = A*P_est*(A.') + V;
end