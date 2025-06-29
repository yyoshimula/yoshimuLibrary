%[text] # precession quaternion from mean to TOD
%[text] `jd0`: Julian day, day, scalar
%[text] `jd1`: Julian day, day, scalar
%[text] `scalar:` specifies the definition of the quaternion 
%[text] `scalar == 0`
%[text] ${\\bf q} = \[q\_0,q\_1,q\_2,q\_3\]^T=\[\\cos(\\theta/2), {\\bf e}^T\\sin(\\theta/2)\]^T$
%[text] `scalar == 4`
%[text] ${\\bf q} = \[q\_1,q\_2,q\_3, q\_4\]^T=\[{\\bf e}^T\\sin(\\theta/2), \\cos(\\theta/2)\]^T$
%[text] const: orbital constants
%[text] q: quaternion from jd0 to jd1 when jd0 = J2000.0 means the rotation from ICRF to TOD frame
%[text] ## note
%[text] NA
%[text] ## references 
%[text] Montenbruck, O., and Gill, E., Satellite Orbits, Berlin, Heidelberg: Springer Science & Business Media, 2012. p.176
%[text] ## revisions
%[text] 20211027  y.yoshimura, y.yoshimula@gmail.com
%[text] See also oe2roe, nutationQ, precession.
function q = precessionQ(jd0, jd1, scalar, const)
% arguments
%     jd0 (:,1) {mustBeNumeric}
%     jd1 (:,1) {mustBeNumeric}
%     scalar (1,1) {mustBeMember(scalar, [0, 4])}
%     const
% end

[zeta, z, theta, ~, ~, ~] = precession(jd0, jd1, const);

the_q = -pi / 2 - z;
N = length(the_q);
q1 = [ zeros(N,2) sin(the_q/2) cos(the_q/2)];

the_q = theta;
q2 = [sin(the_q/2) zeros(N,2) cos(the_q/2)];

the_q = pi / 2 - zeta;
q3 = [zeros(N,2) sin(the_q/2) cos(the_q/2)];

tmp_q = qMult(4, 1, q1, qMult(4, 1, q2, q3));

q = (scalar == 4 ) .* tmp_q...
     + (scalar == 0) .* [tmp_q(:,4), tmp_q(:,1), tmp_q(:,2), tmp_q(:,3)];

end

%[appendix]{"version":"1.0"}
%---
