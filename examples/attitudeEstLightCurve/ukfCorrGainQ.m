%[text] # Correlated covariances and Kalman gain
%[text]  x\_est: state vector: nx1 vector
%[text] X0, sigma points, 1xn vector
%[text] X: sigma points: 2n x n matrix 
%[text] `yEst`, a priori estimated measurement vector: mx1 vector 
%[text] `Y0`, measurement simga points: mx1 vector 
%[text] `Y`, measurement simga points: 2n x m matrix
%[text] `w0c, wic`: weight for covariance, scalarl
%[text] `R`, measurement noise matrix 
%[text] `Pyy`: measurement covarince: mxm matrix 
%[text] `Pxy`, correlated covariance: nxm matrix
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20210209  y.yoshimura, y.yoshimula@gmail.com
function [Pyy, Pxy, K] = ukfCorrGainQ(scalar, xEst, X0, X, yEst, Y0, Y, w0c, wic, R)

yEst = yEst(:); % mx1 vector
Y0 = Y0(:);
m = length(yEst);

Pyy = zeros(m,m);
%[text] ## covariance
Pyy0 = w0c * (Y0 - yEst) * (Y0 - yEst)'; % mxm

for i = 1:6
    Pyy = Pyy + wic .* (Y(i,:)' - yEst) * (Y(i,:)' - yEst)';
end
Pyy = Pyy0 + Pyy + R;

%[text] ## cross correlation
ex = logQ(scalar, qErr(4, X0, xEst)); % 1 x 3
ex = logQ(scalar, qErr(4, xEst,X0)); % 1 x 3

Pxy0 = w0c * (ex' * (Y0 - yEst)');

Pxy = zeros(3, m);
for i = 1:6
    ex = logQ(scalar, qErr(4, X(i,:), xEst));        
    Pxy = Pxy + wic .* (ex' * (Y(i,:)' - yEst)');
end
Pxy = Pxy0 + Pxy;
%[text] ## Kalman gain
K = Pxy * Pyy^(-1);

end


%[appendix]{"version":"1.0"}
%---
