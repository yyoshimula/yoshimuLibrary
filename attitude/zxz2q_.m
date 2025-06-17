%[text] # ZXZ Euler angle to quaternion
%[text] `phi`: rotation angle around Z-axis (1st rotation)
%[text] `theta`: rotation angle around X-axis (2nd rotation)
%[text] `psi`: rotation angle around Z-axis (3rd rotation)
%[text] `scalar:` specifies the definition of the quaternion 
%[text] `scalar == 0`
%[text] ${\\bf q} = \[q\_0,q\_1,q\_2,q\_3\]^T=\[\\cos(\\theta/2), {\\bf e}^T\\sin(\\theta/2)\]^T$
%[text] `scalar == 4`
%[text] ${\\bf q} = \[q\_1,q\_2,q\_3, q\_4\]^T=\[{\\bf e}^T\\sin(\\theta/2), \\cos(\\theta/2)\]^T$
%[text] `q`: quaternions, nx4 matrix
%[text] `p`: quaternions, nx4 matrix 
%[text] `output`: quaternions, nx4 matrix 
%[text] ## note
%[text] NA
%[text] ## references 
%[text] Markley, F. L., and Crassidis, J. L., Fundamentals of Spacecraft Attitude Determination and Control, New York, NY: Springer, 2014. 
%[text] ## revisions
%[text] 20150101  y.yoshimura, y.yoshimula@gmail.com
%[text] See also qMult.
function q = zxz2q_(scalar, phi, theta, psi)

n = size(phi,1);
%[text] ### each quaternion around Z-axis, X-axis, and Z-axis
q1 = [zeros(n,2), sin(phi./2), cos(phi./2)]; % nx4 matrix
q2 = [sin(theta./2), zeros(n,2), cos(theta./2)];
q3 = [zeros(n,2), sin(psi./2), cos(psi./2)];
tmp = qMult(4, 1, q3, qMult(4,1,q2,q1)); % nx4 matrix
%[text] ## quaternion
q = (scalar == 0) .* [tmp(:,4), tmp(:,1:3)] ...
    + (scalar == 4) .* tmp;

end

%[appendix]{"version":"1.0"}
%---
