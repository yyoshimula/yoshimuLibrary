%[text] # covariance using sigma points
%[text] `xEst:` state vector: 1xn vector
%[text] `X`, sigma points: (2n+1) x n matrix
%[text] `wc:` weight for covariance
%[text] `Q:` process noise matrix, n x n matrix
%[text] `Pcov:`, a priori covariance: n x n matrix
%[text] ## note
%[text] NA
%[text] ## references
%[text] NA
%[text] ## revisions
%[text] 20210209  y.yoshimura, y.yoshimula@gmail.com
%[text] See also ukfSigma, ukf.
function Pcov = ukfCov(xEst, X, wc, Q)

xEst = xEst(:); % column vector
n = length(xEst);
Pcov = zeros(n,n);

%[text] ## covariance sum
for i = 1:2*n+1
    Pcov = Pcov + wc(i) * (X(i,:)' - xEst) * (X(i,:)' - xEst)';
end

Pcov = Pcov + Q;

end

%[appendix]{"version":"1.0"}
%---
