%[text] # nutation quaternion from mean of date (MOD) to true of date(TOD)
%[text] ## inputs
%[text] `jd`: Julian day, day, nx1
%[text] `scalar:` specifies the definition of the quaternion 
%[text] `scalar == 0`
%[text] ${\\bf q} = \[q\_0,q\_1,q\_2,q\_3\]^T=\[\\cos(\\theta/2), {\\bf e}^T\\sin(\\theta/2)\]^T$
%[text] `scalar == 4`
%[text] ${\\bf q} = \[q\_1,q\_2,q\_3, q\_4\]^T=\[{\\bf e}^T\\sin(\\theta/2), \\cos(\\theta/2)\]^T$
%[text] const: orbital constants
%[text] ## output
%[text] q: quaternion from mean to TOD, nx4 matrix
%[text] ## note
%[text] NA
%[text] ## references 
%[text] Explanatory Supplement to the Astronomical Almanac, p. 114. Satellite Orbits, p.180
%[text] ## revisions
%[text] 20211027  y.yoshimura, y.yoshimula@gmail.com
%[text] See also oe2roe, orbitConst.
function q = nutationQ(jd, scalar, const)
% arguments
%     jd (:,1) {mustBeNumeric}
%     scalar (1,1) {mustBeMember(scalar, [0, 4])}
%     const
% end

%[text] ## mean obliquity of the ecliptic at the epoch, rad
e = obliquity(jd);

% dPsi: nutation in longitude, rad
% dEpsi: nutation in obliquity, rad
[dPsi, dEpsi] = nutation(jd,const);

ePrime = dEpsi + e;

%[text] R1 = DCM1axis(1, -e\_prime);
%[text] R2 = DCM1axis(3, -dPsi);
%[text] R3 = DCM1axis(1, e);
%[text] dcm = R1 \* R2 \* R3;
%[text] 
%[text] forth part of quaternion is assumed to be scalar part, i.e., q = \[e.\*sin(theta/2); cos(theta/2)\]
N = length(jd);

q1 = [ones(N,1).*sin(-ePrime./2), zeros(N,2), cos(-ePrime./2)];

q2 = [zeros(N,2), ones(N,1).*sin(-dPsi./2), cos(-dPsi./2)];

q3 = [ones(N,1).*sin(e./2), zeros(N,2), cos(e./2)];

tmpQ = qMult(4, 1, q1, qMult(4, 1, q2, q3));

q = (scalar == 4 ) .* tmpQ...
     + (scalar == 0) .* [tmpQ(:,4), tmpQ(:,1), tmpQ(:,2), tmpQ(:,3)];

end

%[appendix]{"version":"1.0"}
%---
