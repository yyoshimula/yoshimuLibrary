function Pcov = ukfCovQ(scalar, xEst, X0, X, w0c, wic, Q)

Pmat = zeros(3,3);

% covariance sum
e = logQ(scalar, qErr(scalar, X0, xEst)); % 1 x 3

P0 = w0c .* (e' * e);

e = logQ(scalar, qErr(scalar, X, repmat(xEst, 6, 1))); % 2n x 3

for i = 1:6
    Pmat = Pmat + wic .* e(i,:)' * e(i,:);
end

Pcov = P0 + Pmat + Q;

end
