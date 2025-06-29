%[text] この関数の簡単な概要です。
%[text] この関数の詳細な説明です。
function qMat = qExp(scalar, dt, w)
arguments (Input)
    scalar
    dt
    w
end

arguments (Output)
    qMat
end

w = w(:);
wNorm = norm(w);
qMat = eye(4,4) .* cos(0.5 * dt * wNorm) + qMultMat(scalar, 1, [w', 0]) ./ wNorm .* sin(0.5 * dt * wNorm);


end

%[appendix]{"version":"1.0"}
%---
