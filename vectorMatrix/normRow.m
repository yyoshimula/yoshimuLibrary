%[text] # Normalize row vectors of matrix

function B = normRow(A)

Anorm = vecnorm(A, 2, 2);

if (Anorm > eps)
    B = A ./ Anorm;
else
    B = A;
end

end

%[appendix]{"version":"1.0"}
%---
