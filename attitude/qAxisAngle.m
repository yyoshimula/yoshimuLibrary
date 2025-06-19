%[text] # Euler axis and rotation angle from quaternions
%[text] `scalar:` specifies the definition of the quaternion 
%[text] `scalar == 0`
%[text] ${\\bf q} = \[q\_0,q\_1,q\_2,q\_3\]^T=\[\\cos(\\theta/2), {\\bf e}^T\\sin(\\theta/2)\]^T$
%[text] `scalar == 4`
%[text] ${\\bf q} = \[q\_1,q\_2,q\_3, q\_4\]^T=\[{\\bf e}^T\\sin(\\theta/2), \\cos(\\theta/2)\]^T$
%[text] `q`: quaternions, nx4 matrix
%[text] `eAxis`: Euler axis, $\\bf e$, nx3 matrix
%[text] `eAngle`: rotation angle, $\\theta \\in \[0, 2\\pi)$, nx1 matrix, rad,
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20150101  y.yoshimura, y.yoshimula@gmail.com
%[text] See also qCon, qInv.
function [eAxis, eAngle] = qAxisAngle(scalar, q)
% arguments
%     scalar (1,1) {mustBeMember(scalar, [0, 4])}
%     q (:,4) {mustBeNumeric}
% end

if scalar == 0
    eAngle = 2.0 .* acos(q(:,1));    
    eAxis = q(:,2:4) ./ sin(eAngle ./ 2.0);
    
elseif scalar == 4
    eAngle = 2.0 .* acos(q(:,4));
    eAxis = q(:,1:3) ./ sin(eAngle ./ 2.0);

else
    error('definition of quaternions is unclear');
    
end

eAngle = mod(eAngle, 2*pi);

end

%[appendix]{"version":"1.0"}
%---
