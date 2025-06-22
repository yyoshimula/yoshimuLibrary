%[text] # calculating sigma points
%[text] lam: tuning parameter of UKF
%[text] P: covariance matrix, nxn matrix
%[text] x, state vector, nx1 or 1xn vector
%[text] method:     `'svd'` (default) singular value decomposition for P sqaure root
%[text]  `'chol'`: Cholesky decomposition for P square root
%[text] `X:` sigma points: (2n+1) x n matrix
%[text] ## note
%[text] NA
%[text] ## references
%[text] NA
%[text] ## revisions
%[text] 20240809 arguments modified
%[text] 20181210  y.yoshimura, y.yoshimula@gmail.com
%[text] See also ukfCov.
function X = ukfSigma(lam, P, x, method)
% method引数が設定されていない場合はデフォルト値を使用
if nargin < 4
    method = 'svd';
end

n = length(x);
%[text] ### square root matrix

if strcmp(method, 'svd') % SVD
    [U, S, ~] = svd(P);
    Psq = U * sqrt(S);
    
elseif strcmp(method, 'chol')
    Psq = chol(P)'; % = sqrt(P)
else
    error('no proper decomposition method, use svd or chol')
end

%[text] ### sigma points
sig = [sqrt(n+lam).*Psq -sqrt(n+lam).*Psq]; % n x 2n matrix
X = [x;
    x + sig']; % (2n+1) x n matrix

end

%[appendix]{"version":"1.0"}
%---
