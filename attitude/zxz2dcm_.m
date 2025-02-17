% ZXZ Euler angle to directional cosine matrix (rotation matrix)
% phi: rotation angle around Z-axis (1st rotation)
% theta: rotation angle around X-axis (2nd rotation)
% psi: rotation angle around Z-axis (3rd rotation)
% R: directional cosine matrix (rotation matrix), 3x3
% note
% NA
% references 
% Markley, F. L., and Crassidis, J. L., Fundamentals of Spacecraft Attitude Determination and Control, New York, NY: Springer, 2014. 
% revisions
% 20150101  y.yoshimura, y.yoshimula@gmail.com
% See also q2DCM, zyx2dcm.

function R = zxz2dcm_(phi, theta, psi)

R = dcm1axisZ_(psi) * dcm1axisX_(theta) * dcm1axisZ_(phi);

end


