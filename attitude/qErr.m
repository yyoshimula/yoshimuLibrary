%[text] # calculate quaternion error 
%[text] $q\_d$を回転させて$q$に一致させるためのerror quaternion
%[text] `scalar:` specifies the definition of the quaternion 
%[text] `scalar == 0`
%[text] ${\\bf q} = \[q\_0,q\_1,q\_2,q\_3\]^T=\[\\cos(\\theta/2), {\\bf e}^T\\sin(\\theta/2)\]^T$
%[text] `scalar == 4`
%[text] ${\\bf q} = \[q\_1,q\_2,q\_3, q\_4\]^T=\[{\\bf e}^T\\sin(\\theta/2), \\cos(\\theta/2)\]^T$
%[text] `q`: quaternions, nx4 matrix
%[text] qd: desired (or estimated ) quaternion, nx4 matrix 
%[text] `output`: quaternion error, nx4 matrix 
%[text] ${\\bf q}\_e = {\\bf q}\_d^{-1} \\odot {\\bf q} = {\\bf q} \\otimes {\\bf q}\_d^{-1} $
%[text] ## note
%[text] NA
%[text] ## references 
%[text] Markley, F. L., and Crassidis, J. L., Fundamentals of Spacecraft Attitude Determination and Control, New York, NY: Springer, 2014. 
%[text] ## revisions
%[text] 20150101  y.yoshimura, y.yoshimula@gmail.com
%[text] See also qCon, qInv.
function qe = qErr(scalar, q, qd)

% q \otimes p = [ q0 * p0 - qv' * pv
%               q0 .* pv + p0 .* qv - cross(qv,pv)]

qInvD = qInv(scalar, qd);

qe = qMult(scalar, 1, q, qInvD);

end

%[appendix]{"version":"1.0"}
%---
