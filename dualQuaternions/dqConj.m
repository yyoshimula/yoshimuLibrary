%[text] # Dual quaternion conjugate
%[text] 共役dua lquaternion
%[text] `scalar:` specify the definition of the quaternion 
%[text] `scalar == 0`
%[text] ${\\bf q} = \[q\_0,q\_1,q\_2,q\_3\]^T=\[\\cos(\\theta/2), {\\bf e}^T\\sin(\\theta/2)\]^T$
%[text] `scalar == 4`
%[text] ${\\bf q} = \[q\_1,q\_2,q\_3, q\_4\]^T=\[{\\bf e}^T\\sin(\\theta/2), \\cos(\\theta/2)\]^T$
%[text] dq: dual quaternion \[real, dual\], nx8
%[text] dqCon: dual quaternion conjugate \[real, dual\]，nx8
%[text] ## note
%[text] NA
%[text] ## references 
%[text] Sveier, A., & Egeland, O. (2020). Dual Quaternion Particle Filtering for Pose Estimation. IEEE Transactions on Control Systems Technology, 1-14.
%[text] ## revisions
%[text] 20220526  y.yoshimura
%[text] See also qMult.
function invDq = dqConj(scalar, dq)
% arguments
%     scalar (1,1) {mustBeMember(scalar, [0, 4])}
%     dq (:,8) {mustBeNumeric}
% end

if scalar == 0
    invDq = [dq(:,1) -dq(:,2:4) dq(:,5) -dq(:,6:8)];
else
    invDq = [-dq(:,1:3) dq(:,4) -dq(:,5:7) dq(:,8)];
end

end


%[appendix]{"version":"1.0"}
%---
