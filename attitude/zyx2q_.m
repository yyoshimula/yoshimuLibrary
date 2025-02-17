% ZYX Euler angle to quaternion
% phi: rotation angle around Z-axis (1st rotation)
% theta: rotation angle around Y-axis (2nd rotation)
% psi: rotation angle around X-axis (3rd rotation)
% scalar: specifies the definition of the quaternion 
% scalar == 0
% 
% scalar == 4
% 
% q: quaternions, nx4 matrix
% p: quaternions, nx4 matrix 
% output: quaternions, nx4 matrix 
% note
% NA
% references 
% Markley, F. L., and Crassidis, J. L., Fundamentals of Spacecraft Attitude Determination and Control, New York, NY: Springer, 2014. 
% revisions
% 20150101  y.yoshimura, y.yoshimula@gmail.com
% See also qMult.

function q = zyx2q_(scalar, phi, theta, psi)

n = size(phi,1);

% each quaternion around Z-axis, Y-axis, and X-axis
% q4: scalar partで計算

qPhi = [zeros(n,2), sin(phi/2), cos(phi/2)];
qTheta = [zeros(n,1), sin(theta/2), zeros(n,1), cos(theta/2)];
qPsi = [sin(psi/2), zeros(n,2), cos(psi/2)];
q = qMult(4,1, qPsi, qMult(4,1, qTheta, qPhi));

% quaternion for output
q = (scalar == 0) .* [q(:,4), q(:,1:3)] ...
    + (scalar == 4) .* q;

end
