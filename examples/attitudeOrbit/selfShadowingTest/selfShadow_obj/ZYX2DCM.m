function R = ZYX2DCM(phi, theta, psi)
%  2012/2/20 yy  ZYX Euler angles to rotation matrix
%  input : ZYX euler angles(3x1 vector), [phi, theta, psi] [rad]
%  output: Rotation matrix from inertial coordinate to body-fixed one
% ----------------------------------------------------
% DCM denotes a rotation from inertial frame to body frame
% i.e. b = R*i
% ZYX Euler angles
% Euler[phi (around z-axis); theta (around y-axis); psi (around x-axis)]
% R = ZYX2DCM(phi, theta, psi)

R = [cos(theta)*cos(phi) cos(theta)*sin(phi) -sin(theta)
    sin(theta)*cos(phi)*sin(psi)-sin(phi)*cos(psi) sin(theta)*sin(phi)*sin(psi)+cos(phi)*cos(psi) cos(theta)*sin(psi)
    sin(theta)*cos(phi)*cos(psi)+sin(phi)*sin(psi) sin(theta)*sin(phi)*cos(psi)-cos(phi)*sin(psi) cos(theta)*cos(psi)];

end