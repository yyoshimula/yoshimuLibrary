%[text] # calculating ZYZ (3-2-3) Euler angles from quaternions 
%[text] `scalar:` specify the definition of the quaternion 
%[text] `scalar == 0`
%[text] ${\\bf q} = \[q\_0,q\_1,q\_2,q\_3\]^T=\[\\cos(\\theta/2), {\\bf e}^T\\sin(\\theta/2)\]^T$
%[text] `scalar == 4`
%[text] ${\\bf q} = \[q\_1,q\_2,q\_3, q\_4\]^T=\[{\\bf e}^T\\sin(\\theta/2), \\cos(\\theta/2)\]^T$
%[text] `q`: quaternions, nx4 matrix
%[text] `euler:` Euler angles, \[$\\phi,\\theta,\\psi$\], rad, nx3 matrix
%[text] where  $\\phi\n$: 1st rotation around z-axis
%[text]  $\\theta\n$: 2nd rotation around y-axis 
%[text]  $\\psi$: 3rd rotation around z-axis 
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20230809  y.yoshimura, y.yoshimula@gmail.com
%[text] See also q2zyx.
function euler = q2zyz(scalar, q)
% arguments
%     scalar (1,1) {mustBeMember(scalar, [0, 4])}
%     q (:,4) {mustBeNumeric}
% end

% calculate DCM using the definition that q4 = cos(theta/2)
if scalar == 0
    q = [q(:,2), q(:,3),q(:,4), q(:,1)];
end

R23 = 2.0 * (q(:,2) .* q(:,3) + q(:,4) .* q(:,1));
R13 = 2.0 * (q(:,1) .* q(:,3) - q(:,4) .* q(:,2));

R31 = 2.0 * (q(:,1) .* q(:,3) + q(:,4) .* q(:,2));
R32 = 2.0 * (q(:,2) .* q(:,3) - q(:,4) .* q(:,1));

R33 = -q(:,1).^2 - q(:,2).^2 + q(:,3).^2 + q(:,4).^2;

psi   = atan2(R23, -R13);
phi   = atan2(R32, R31);
theta = atan2(sqrt(R31.^2 + R32.^2), R33);

euler = [phi, theta, psi]; % [rad]

end


%[appendix]{"version":"1.0"}
%---
