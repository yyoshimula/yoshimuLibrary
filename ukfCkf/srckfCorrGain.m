%[text] # Correlated covariances and Kalman gain for (spherical radial) cubature Kalman filter
%[text] ## input
%[text] xEst: state vector: nx1 vector
%[text] X: sigma cubature points: 2n x n matrix 
%[text] `yEst`, a priori estimated measurement vector: mx1 vector 
%[text] `Y`, measurement simga cubature points: 2n x m matrix
%[text] `wc`, weight for covariance: 2n+1 x 1 vector
%[text] `R`, measurement noise matrix 
%[text] ## output
%[text] `Pyy`: measurement covarince: mxm matrix 
%[text] `Pxy`, correlated covariance: nxm matrix
%[text] K: Kalman gain
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20210209  y.yoshimura, y.yoshimula@gmail.com
function [Pyy, Pxy, K] = srckfCorrGain(xEst, X, yEst, Y, wc, R)

nSigma = size(X, 1);
n = length(xEst);
m = length(yEst);

Pyy = zeros(m, m);
Pxy = zeros(n, m);

for i = 1:nSigma
    dy = Y(i,:)' - yEst(:);
    dx = X(i,:)' - xEst(:);
    Pyy = Pyy + wc(i) * (dy * dy');
    Pxy = Pxy + wc(i) * (dx * dy');
end

Pyy = Pyy + R;
K = Pxy * Pyy^(-1);


end



%[appendix]{"version":"1.0"}
%---
