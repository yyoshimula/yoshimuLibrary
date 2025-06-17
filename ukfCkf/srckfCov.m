%[text] # covariance using (spherical radial) sigma cubature points
%[text] ## input
%[text] `xEst:` state vector: 1xn vector
%[text] `X`, sigma cubature points: 2n x n matrix 
%[text] `Q:` process noise matrix, n x n matrix
%[text] ## output
%[text] `Pout:`, a priori covariance: n x n matrix
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20210209  y.yoshimura, y.yoshimula@gmail.com
%[text] See also ukfSigma, ukf.
function Pout = srckfCov(xEst, X, wc, Q)

nSigma = size(X, 1);
n = length(xEst);
Pout = zeros(n, n);

for i = 1:nSigma
    dx = X(i,:)' - xEst(:);
    Pout = Pout + wc(i) * (dx * dx');
end

Pout = Pout + Q;

end



%[appendix]{"version":"1.0"}
%---
