%[text] # decoupled translational and rotational motion
%[text] `scalar:` specify the definition of the quaternion 
%[text] `scalar == 0`
%[text] ${\\bf q} = \[q\_0,q\_1,q\_2,q\_3\]^T=\[\\cos(\\theta/2), {\\bf e}^T\\sin(\\theta/2)\]^T$
%[text] `scalar == 4`
%[text] ${\\bf q} = \[q\_1,q\_2,q\_3, q\_4\]^T=\[{\\bf e}^T\\sin(\\theta/2), \\cos(\\theta/2)\]^T$
%[text] x: dual quaternion \[real, dual\], 1x8
%[text] ## note
%[text] NA
%[text] ## references 
%[text] Sveier, A., & Egeland, O. (2020). Dual Quaternion Particle Filtering for Pose Estimation. IEEE Transactions on Control Systems Technology, 1-14.
%[text] ## revisions
%[text] 20231219  y.yoshimura y.yoshimula@gmail.com
%[text] See also pos2dq.
function dxdt = decEOM(~, x_)
global m MOI mu

q = x_(1:4);
w = x_(5:7);
r = x_(8:10);
v = x_(11:13);

a = -mu / norm(r)^3 .* r;

dxdt = [qKine(4, q', w')'
    MOI^(-1) * (-cross(w,MOI*w))
    v
    a];


end

%[appendix]{"version":"1.0"}
%---
