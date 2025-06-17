%[text] # calculating sigma cubature points
%[text] ## input
%[text] P: covariance matrix, nxn matrix
%[text] x, state vector, nx1 or 1xn vector 
%[text] ## output
%[text] `Xout:` cubature sigma points: 2n x n matrix 
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20250617  y.yoshimura, y.yoshimula@gmail.com
%[text] See also ckfCov.
function Xout = srckfSigma(P, x)

x = x(:); % column vector
n = length(x);

% Cholesky分解
S = chol(P)';

% 球面シグマ点
xi = sqrt(n) * S;  % n x n matrix

% シグマ点構成
Xout = [x, x + xi, x - xi]; % n x (2n+1)
Xout = Xout';  % (2n+1) x n

end



%[appendix]{"version":"1.0"}
%---
