%[text] # Correlated covariances and Kalman gain
%[text]  x\_est: state vector: nx1 vector
%[text] X: sigma points: 2n+1 x n matrix 
%[text] `yEst`, a priori estimated measurement vector: mx1 vector 
%[text] `Y`, measurement simga points: 2n+1 x m matrix
%[text] `wc`: weight for covariance, scalar or 1x2n+1 vector
%[text] `R`, measurement noise matrix 
%[text] `Pyy`: measurement covarince: mxm matrix 
%[text] `Pxy`, correlated covariance: nxm matrix
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20210209  y.yoshimura, y.yoshimula@gmail.com
function [Pyy, Pxy, K] = ukfCorrGain(xEst, X, yEst, Y, wc, R)

xEst = xEst(:); % column vector
n = length(xEst);

yEst = yEst(:); % mx1 vector
m = length(yEst);

Pyy = zeros(m,m);
%[text] ## covariance

for i = 1:(2*n+1)
    Pyy = Pyy + wc(i) .* (Y(i,:)' - yEst) * (Y(i,:)' - yEst)';
end
Pyy = Pyy + R;

%[text] ## cross correlation
Pxy = zeros(n, m);
for i = 1:(2*n+1)
    Pxy = Pxy + wc(i) .* (X(i,:)' - xEst) * (Y(i,:)' - yEst)';
end

%[text] ## Kalman gain
K = Pxy * Pyy^(-1);

end

%[appendix]{"version":"1.0"}
%---
