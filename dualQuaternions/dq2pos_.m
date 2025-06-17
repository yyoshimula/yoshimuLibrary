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
function [r, q] = dq2pos_(inertial, scalar, dq)

%[text] ### quaternion = real part
q = dq(:,1:4);

%[text] ### dual part
tmpR = (inertial == 0) .* qMult_(scalar, 1, dq(:,5:8), qInv_(scalar,q)) ...
    + (inertial == 1) .* qMult_(scalar, 1, qInv_(scalar, q), dq(:,5:8));
tmpR = tmpR .* 2.0;

r = (scalar == 0) .* tmpR(:,2:4) + (scalar == 4) .* tmpR(:,1:3);

end


%[appendix]{"version":"1.0"}
%---
