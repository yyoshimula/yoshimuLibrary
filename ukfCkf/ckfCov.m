%[text] # covariance using sigma cubature points
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
function Pout = ckfCov(xEst, X, Q)

n = size(X,2);
P = zeros(n,n);

%[text] ## covariance sum
for i = 1:2*n
    P = P + X(i,:)' * X(i,:);
end
P = P ./ (2 * n);

Pout = P - xEst' * xEst + Q;

end


%[appendix]{"version":"1.0"}
%---
