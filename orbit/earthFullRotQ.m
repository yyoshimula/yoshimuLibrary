%[text] # Earth rotation including precession and nutation
%[text] `jd0`: Julian day, day, scalar
%[text] `jd1`: Julian day, day, nx1 vector
%[text] `scalar:` specifies the definition of the quaternion 
%[text] `scalar == 0`
%[text] ${\\bf q} = \[q\_0,q\_1,q\_2,q\_3\]^T=\[\\cos(\\theta/2), {\\bf e}^T\\sin(\\theta/2)\]^T$
%[text] `scalar == 4`
%[text] ${\\bf q} = \[q\_1,q\_2,q\_3, q\_4\]^T=\[{\\bf e}^T\\sin(\\theta/2), \\cos(\\theta/2)\]^T$
%[text] q: quaternion from epoch jd0 to jd1, nx4 matrix
%[text] ## note
%[text] Polar motion is not included.
%[text] ## references
%[text] NA
%[text] ## revisions
%[text] 20211027  y.yoshimura, y.yoshimula@gmail.com
%[text] See also nutationQ, precessionQ.
function q = earthFullRotQ(jd0, jd1, scalar, const)
arguments
    jd0 (:,1) {mustBeNumeric}
    jd1 (:,1) {mustBeNumeric}
    scalar (1,1) {mustBeMember(scalar, [0, 4])}
    const
end

%[text] ## Precession 
qP = precessionQ(jd0, jd1, scalar, const);

%[text] ## Nutation
qN = nutationQ(jd1, scalar, const);

%[text] ## TOD = Nutation \* Precession \* J2000
q = qMult(scalar, 1, qN, qP);

end

%[appendix]{"version":"1.0"}
%---
