function [Pest, xEst] = ukf(f, h, B, Q, R, lam, P, x, y, wm, wc, wic)
% ----------------------------------------------------------------------
%   Unscented Kalman Filter
%    20181210  y.yoshimura
%    Inputs: f, system descrete dynamics: x(k+1) = f(x(k)), nx1 vector
%            h, observation: y(k) = h(x(k)), mx1 vector
%            B, system noise sensitivity matrix, n x r matrix
%            Q, system noise: rx1 vector
%            R, observation noise: mx1 vector
%            lam, tuning parameter for UKF
%            P, covariance matrix, nxn matrix
%            x, state vector, nx1 vector
%            y, observation vector, mx1 vector
%   
%   Outputs: P_est, estimated covariance matrix, nxn matrix
%            x_est, estimated state, nx1 vector
%   related function files: ukf_sigma, ukf_cov, ukf_corr, testUKF
%   note: 
%   cf:
%   revisions;
%
%   (c) 2018 yasuhiro yoshimura
%----------------------------------------------------------------------
% sigma points
X = ukfSigma(lam, P, x);

% sigma point propagation
X = f(X');

X = X';
% a priroi estimation, 1xn vector
xEst = sum(wim .* X, 1);

% covariance
P = ukfCov(xEst, X, wic, Q);

% output sigma points
Y = h(X);
% a priori output estimation, mx1 vector
yEst = sum(wm .* Y, 1);

% Calculate correlation
[Pyy, ~, K] = ukfCorrGain(xEst, X, yEst, Y, wc, R);

% Gain and Update
Pest = P - K * Pyy * K';
xEst = xEst + (K * (y - yEst))';

end