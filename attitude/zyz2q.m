%[text] # ZYZ (3-2-3) Euler angle to quaternion
%[text] ## input
%[text] `phi`: rotation angle around Z-axis (1st rotation)
%[text] `theta`: rotation angle around Y-axis (2nd rotation)
%[text] `psi`: rotation angle around Z-axis (3rd rotation)
%[text] `scalar:` specifies the definition of the quaternion 
%[text] `scalar == 0`
%[text] ${\\bf q} = \[q\_0,q\_1,q\_2,q\_3\]^T=\[\\cos(\\theta/2), {\\bf e}^T\\sin(\\theta/2)\]^T$
%[text] `scalar == 4`
%[text] ${\\bf q} = \[q\_1,q\_2,q\_3, q\_4\]^T=\[{\\bf e}^T\\sin(\\theta/2), \\cos(\\theta/2)\]^T$
%[text] `q`: quaternions, nx4 matrix
%[text] `p`: quaternions, nx4 matrix 
%[text] ## output
%[text] q: quaternions, nx4 matrix 
%[text] ## note
%[text] NA
%[text] ## references 
%[text] Markley, F. L., and Crassidis, J. L., Fundamentals of Spacecraft Attitude Determination and Control, New York, NY: Springer, 2014. 
%[text] ## revisions
%[text] 20240315  y.yoshimura, y.yoshimula@gmail.com
%[text] See also qMult, zyx2q.
function q = zyz2q(scalar, phi, theta, psi)
% arguments
%     scalar (1,1) {mustBeMember(scalar, [0, 4])}
%     phi (:,1) {mustBeNumeric}
%     theta (:,1) {mustBeNumeric}
%     psi (:,1) {mustBeNumeric}
% end

n = size(phi,1);
%[text] ### each quaternion around Z-axis, Y-axis, and Z-axis
%[text] q4: scalar partで計算
qPhi = [zeros(n,2), sin(phi/2), cos(phi/2)];
qTheta = [zeros(n,1), sin(theta/2), zeros(n,1), cos(theta/2)];
qPsi = [zeros(n,2), sin(psi/2), cos(psi/2)];
q = qMult(4,1, qPsi, qMult(4,1, qTheta, qPhi));
%[text] ## quaternion
if scalar == 0
    q = [q(:,4), q(:,1:3)];
    
elseif scalar == 4
    % do nothing

else
    error('quaternion definition is unclear')

end

end


%[appendix]{"version":"1.0"}
%---
