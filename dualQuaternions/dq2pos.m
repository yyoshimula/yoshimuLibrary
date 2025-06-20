%[text] # transforming dual quaternion to position and quaternion
%[text] dual quaternionを位置ベクトルとquaternionへ変換
%[text] ## input
%[text] `inertial`: 1: position vector is expressed with inertial frame, 0: body-fixed frame
%[text] `scalar:` specify the definition of the quaternion 
%[text] `scalar == 0`
%[text] ${\\bf q} = \[q\_0,q\_1,q\_2,q\_3\]^T=\[\\cos(\\theta/2), {\\bf e}^T\\sin(\\theta/2)\]^T$
%[text] `scalar == 4`
%[text] ${\\bf q} = \[q\_1,q\_2,q\_3, q\_4\]^T=\[{\\bf e}^T\\sin(\\theta/2), \\cos(\\theta/2)\]^T$
%[text] dq: dual quaternion \[real, dual\], nx8
%[text] ## output
%[text] r: position vector, nx3
%[text] q: quaternion, nx4
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20231219  y.yoshimura
%[text] See also dqMult, pos2dq.
function [r, q] = dq2pos(inertial, scalar, dq)
% arguments
%     inertial (1,1) {mustBeMember(inertial, [0, 1])}
%     scalar (1,1) {mustBeMember(scalar, [0, 4])}
%     dq (:,8) {mustBeNumeric}
% end

%[text] ### quaternion = real part
q = dq(:,1:4);

%[text] ### dual part
if inertial == 0
    tmpR = qMult(scalar, 1, dq(:,5:8), qInv(scalar,q));
elseif inertial == 1
    tmpR = qMult(scalar, 1, qInv(scalar, q), dq(:,5:8));
else
    error('inertial definition error')
end

tmpR = tmpR .* 2.0;

if scalar == 0
    r = tmpR(:,2:4);
elseif scalar == 4
    r = tmpR(:,1:3);
else
    error('quaternion definition error')
end

end

%[appendix]{"version":"1.0"}
%---
