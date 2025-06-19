%[text] # quaternion kinematics
%[text] kinematic equation using quaternions:
%[text] ${\\bf \\dot{q}} = \\frac{1}{2}{\\bf \\omega}\\otimes {\\bf q}$
%[text] `scalar:` specify the definition of the quaternion 
%[text] `scalar == 0`
%[text] ${\\bf q} = \[q\_0,q\_1,q\_2,q\_3\]^T=\[\\cos(\\theta/2), {\\bf e}^T\\sin(\\theta/2)\]^T$
%[text] `scalar == 4`
%[text] ${\\bf q} = \[q\_1,q\_2,q\_3, q\_4\]^T=\[{\\bf e}^T\\sin(\\theta/2), \\cos(\\theta/2)\]^T$
%[text] R, rotation matrix, 3x3 matrix 
%[text] `q`: quaternions, 1x4 vector 
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20200901  y.yoshimura, y.yoshimula@gmail.com
%[text] See also q2DCM.
function qKine = qKine(scalar, q, w)
% arguments 
%     scalar (1,1) {mustBeMember(scalar, [0, 4])}
%     q (:,4) {mustBeNumeric}
%     w (:,3) {mustBeNumeric}
% end

q = q(:);
wx = w(1);
wy = w(2);
wz = w(3);

if (scalar == 0)
    mat = [0 -wx -wy -wz
        wx 0 wz -wy
        wy -wz 0 wx
        wz wy -wx 0];
    qKine = 0.5 .* mat * q;
    
    qKine = qKine';
    
elseif (scalar == 4)
    mat = [0 wz -wy wx
        -wz 0  wx wy
        wy -wx 0 wz
        -wx -wy -wz 0];
    qKine = 0.5 .* mat * q;
    
    qKine = qKine';
    
else
    error('definition of quaternions is unclear');
    
end


end

%[appendix]{"version":"1.0"}
%---
