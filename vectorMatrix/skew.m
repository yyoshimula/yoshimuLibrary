%[text] #  make skew-symmetric matrix from vector x
function S = skew(x)

S = [0 -x(3) x(2)
    x(3) 0 -x(1)
    -x(2) x(1) 0];

end

%[appendix]{"version":"1.0"}
%---
