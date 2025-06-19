%[text] # calculating ZYX Euler angles from quaternions 
%[text] `scalar:` specify the definition of the quaternion 
%[text] `scalar == 0`
%[text] ${\\bf q} = \[q\_0,q\_1,q\_2,q\_3\]^T=\[\\cos(\\theta/2), {\\bf e}^T\\sin(\\theta/2)\]^T$
%[text] `scalar == 4`
%[text] ${\\bf q} = \[q\_1,q\_2,q\_3, q\_4\]^T=\[{\\bf e}^T\\sin(\\theta/2), \\cos(\\theta/2)\]^T$
%[text] `q`: quaternions, nx4 matrix
%[text] `euler:` Euler angles, \[$\\phi,\\theta,\\psi$\], rad, nx3 matrix
%[text] where  $\\phi\n$: 1st rotation around z-axis
%[text]  $\\theta\n$: 2nd rotation around y-axis 
%[text]  $\\psi$: 3rd rotation around x-axis 
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20150101  y.yoshimura, y.yoshimula@gmail.com
%[text] See also q2zyx.
function euler = q2zyx(scalar, q)
% arguments
%     scalar (1,1) {mustBeMember(scalar, [0, 4])}
%     q (:,4) {mustBeNumeric}
% end

% calculate DCM using the definition that q0 = cos(theta/2)
if scalar == 4
    q = [q(:,4), q(:,1), q(:,2), q(:,3)];
end

R23 = 2.0 .* (q(:,3) .* q(:,4) + q(:,1) .* q(:,2));
R33 = q(:,1).^2 - q(:,2).^2 - q(:,3).^2 + q(:,4).^2;
R12 = 2.0 .* (q(:,2) .* q(:,3) + q(:,1) .* q(:,4));
R11 = q(:,1).^2 + q(:,2).^2 - q(:,3).^2 - q(:,4).^2;
R13 = 2.0 .* (q(:,2) .* q(:,4) - q(:,1) .* q(:,3));

psi   = atan2(R23, R33);
phi   = atan2(R12, R11);
theta = atan2(-R13, sqrt(R23.^2 + R33.^2));

euler = [phi, theta, psi]; % [rad]

end

%[appendix]{"version":"1.0"}
%---
