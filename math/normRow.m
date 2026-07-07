%[text] # Normalize row vectors of matrix

function B = normRow(A)

Anorm = vecnorm(A, 2, 2);

% 行ごとに正規化する（ゼロ行はそのまま残す）。旧実装は if (Anorm > eps) が
% 複数行入力で列ベクトル比較となり、ゼロ行が混在すると全行が未正規化のまま返っていた。
B = A;
nz = Anorm > eps;
B(nz, :) = A(nz, :) ./ Anorm(nz);

end

%[appendix]{"version":"1.0"}
%---
