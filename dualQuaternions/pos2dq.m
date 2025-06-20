%[text] # transforming position and quaternion to dual quaternion
%[text] 位置ベクトルとquaternionをdual quaternionへ変換
%[text] `inertial`: 0: body-fixed frame, 1: position vector is expressed with inertial frame
%[text] `scalar:` specify the definition of the quaternion 
%[text] `scalar == 0`
%[text] ${\\bf q} = \[q\_0,q\_1,q\_2,q\_3\]^T=\[\\cos(\\theta/2), {\\bf e}^T\\sin(\\theta/2)\]^T$
%[text] `scalar == 4`
%[text] ${\\bf q} = \[q\_1,q\_2,q\_3, q\_4\]^T=\[{\\bf e}^T\\sin(\\theta/2), \\cos(\\theta/2)\]^T$
%[text] dq: dual quaternion \[real, dual\], nx8
%[text] ## note
%[text] NA
%[text] ## references 
%[text] Sveier, A., & Egeland, O. (2020). Dual Quaternion Particle Filtering for Pose Estimation. IEEE Transactions on Control Systems Technology, 1-14.
%[text] ## revisions
%[text] 20220526  y.yoshimura
%[text] See also qMult.
function dq = pos2dq(inertial, scalar, r, q)
% arguments
%     inertial (1,1) {mustBeMember(inertial, [0, 1])}
%     scalar (1,1) {mustBeMember(scalar, [0, 4])}
%     r (:,3) {mustBeNumeric}
%     q (:,4) {mustBeNumeric}
% end

%[text] ### real part
dqReal = q;

%[text] ### dual part
n = size(r,1);

if scalar == 0
    pos = [zeros(n, 1), r];
else
    pos = [r, zeros(n, 1)];
end

if inertial == 0
    dqDual = 0.5 .* qMult(scalar, 1, pos, q);
else
    dqDual = 0.5 .* qMult(scalar, 1, q, pos);
end

dq = [dqReal, dqDual];


end

%[appendix]{"version":"1.0"}
%---
