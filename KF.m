% Kalman Filter
% Input : A, B, C - System Matrices
% Input : V, W - Process and Measurement Noise Covariance
% Input : x_bar_est : state estimate, y - observation vector
% Input : u - previous input
% Output : x_upd - State Estimate
% Output : P_upd - Posterior State Covariance Matrix
function [x_upd, P_upd] = KF(A,B,C,V,W,P,x_bar_est,y,u)
    n = length(x_bar_est);
    x_prd = A*x_bar_est + B*u; % State Prediction
    P_prd = A*P*(A.') + V; % Prior Covariance
    K = P_prd*(C.')/(C*P_prd*(C.')+W); % Kalman Gain
    x_upd = x_prd + K*(y-C*x_prd); % State Estimate
    P_upd = (eye(n)-K*C)*P_prd; % Posterior Covariance
end
