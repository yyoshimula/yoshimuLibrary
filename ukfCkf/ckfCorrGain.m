%[text] # Correlated covariances and Kalman gain for cubature Kalman filter
%[text] ## input
%[text] xEst: state vector: nx1 vector
%[text] X: sigma cubature points: 2n x n matrix 
%[text] `yEst`, a priori estimated measurement vector: mx1 vector 
%[text] `Y`, measurement simga cubature points: 2n x m matrix
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
function [Pyy, Pxy, K] = ckfCorrGain(xEst, X, yEst, Y, R)
xEst = xEst(:); % column vector
n = length(xEst);

yEst = yEst(:); % mx1 vector
m = length(yEst);

Pyy = zeros(m,m);
%[text] ## covariance
for i = 1:2*n
    Pyy = Pyy + Y(i,:)' * Y(i,:);
end
Pyy = Pyy ./ (2*n) - yEst * yEst' + R;

%[text] ## cross correlation
Pxy = zeros(n, m);
for i = 1:2*n
    Pxy = Pxy + X(i,:)' * Y(i,:);
end
Pxy = Pxy ./ (2*n) - xEst * yEst';
%[text] ## Kalman gain
K = Pxy * Pyy^(-1);

end


%[appendix]{"version":"1.0"}
%---
