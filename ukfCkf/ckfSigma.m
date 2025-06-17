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
function Xout = ckfSigma(P, x)

x = x(:); % column vector
n = length(x);
%[text] ### square root matrix
S = chol(P)'; % i.e., P = S * S'
%[text] ### cubature sigma points
xi = [sqrt(n).*S, -sqrt(n).*S]; % n x 2n matrix
Xout = x + xi; % n x 2n 
Xout = Xout'; % 2n x n matrix

end


%[appendix]{"version":"1.0"}
%---
