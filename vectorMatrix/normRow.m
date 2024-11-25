function B = normRow(A, small)
arguments
    A
    small = 1e-18;
end
% ----------------------------------------------------------------------
%   Normalize row vectors of matrix A
%    20190218  y.yoshimura
%    Inputs: A: matrix, mxn matrix
%   Outputs: B: matrix has unit row vectors, mxn matrix
%   related function files:
%   note:
%   cf:
%   revisions;
%
%   (c) 2019 yasuhiro yoshimura
%----------------------------------------------------------------------
Anorm = vecnorm(A, 2, 2);

if (Anorm > small)
    B = A ./ Anorm;
else
    B = A;
end

end