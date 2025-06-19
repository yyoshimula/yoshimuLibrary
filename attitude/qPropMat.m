%[text] # discrete kinematics of quaternions
%[text] `scalar:` specify the definition of the quaternion 
%[text] `scalar == 0`
%[text] ${\\bf q} = \[q\_0,q\_1,q\_2,q\_3\]^T=\[\\cos(\\theta/2), {\\bf e}^T\\sin(\\theta/2)\]^T$
%[text] `scalar == 4`
%[text] ${\\bf q} = \[q\_1,q\_2,q\_3, q\_4\]^T=\[{\\bf e}^T\\sin(\\theta/2), \\cos(\\theta/2)\]^T$
%[text] `dt`: time span, scalar, s
%[text] `w`: angular rate, 1x3 vector
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20150101  y.yoshimura, y.yoshimula@gmail.com
%[text] See also qInv.
function qProp = qPropMat(scalar, dt, w)
% arguments
%     scalar (1,1) {mustBeMember(scalar, [0, 4])}
%     dt (1,1) {mustBeNumeric}
%     w (1,3) {mustBeNumeric}
% end

w = w(:);

psi = sin(0.5 * norm(w) * dt) .* w ./ norm(w);

if scalar == 0
    qProp = [cos(0.5 * norm(w) * dt), -psi'
        psi, cos(0.5*norm(w)*dt).*eye(3)-skew(psi)];
elseif scalar == 4
    qProp = [cos(0.5*norm(w)*dt).*eye(3)-skew(psi), psi
        -psi', cos(0.5 * norm(w) * dt)];
else
    error('definition of quaternions is unclear');
end

end

%[appendix]{"version":"1.0"}
%---
