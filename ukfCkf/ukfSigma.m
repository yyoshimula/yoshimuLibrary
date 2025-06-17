%[text] # calculating sigma points
%[text] lam: tuning parameter of UKF 
%[text] P: covariance matrix, nxn matrix
%[text] x, state vector, nx1 or 1xn vector 
%[text] method:     `'svd'` (default) singular value decomposition for P sqaure root 
%[text]  `'chol'`: Cholesky decomposition for P square root
%[text] `X0`, sigma points: 1xn vector 
%[text] `X:` sigma points: 2n x n matrix 
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20240809 arguments modified
%[text] 20181210  y.yoshimura, y.yoshimula@gmail.com
%[text] See also ukfCov.
function [X0, X] = ukfSigma(lam, P, x, method)
arguments
    lam
    P
    x
    method = 'svd'
end

x = x(:); % column vector
n = length(x);
%[text] ### square root matrix
if strcmp(method, 'svd') % SVD
    [U, S, ~] = svd(P);
    P_sq = U * sqrt(S);
elseif strcmp(matSqrt, 'chol')
    P_sq = chol(P)'; % = sqrt(P)
else
    error('set proper decomposition method, svd or chol')
end

%[text] ### sigma points
sig = [sqrt(n+lam).*P_sq -sqrt(n+lam).*P_sq]; % n x 2n matrix

X0 = x'; % 1xn vector
X = X0 + sig'; % 2n x n matrix

end

%[appendix]{"version":"1.0"}
%---
