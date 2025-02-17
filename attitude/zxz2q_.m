% ZXZ Euler angle to quaternion
% phi: rotation angle around Z-axis (1st rotation)
% theta: rotation angle around X-axis (2nd rotation)
% psi: rotation angle around Z-axis (3rd rotation)
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

function q = zxz2q_(scalar, phi, theta, psi)

n = size(phi,1);

% each quaternion around Z-axis, X-axis, and Z-axis

q1 = [zeros(n,2), sin(phi./2), cos(phi./2)]; % nx4 matrix
q2 = [sin(theta./2), zeros(n,2), cos(theta./2)];
q3 = [zeros(n,2), sin(psi./2), cos(psi./2)];
tmp = qMult_(4, 1, q3, qMult_(4,1,q2,q1)); % nx4 matrix

% quaternion
q = (scalar == 0) .* [tmp(:,4), tmp(:,1:3)] ...
    + (scalar == 4) .* tmp;

end
