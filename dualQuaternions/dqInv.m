%[text] # Inverse dual quaternions
%[text] 逆dual quaternion
%[text] scalar: 0: quaternionのスカラ部q0としてq=\[q0 q1 q2 q3\]という定義
%[text]             4: quaternionのスカラ部q4としてq=\[q1 q2 q3 q4\]という定義
%[text] dq: dual quaternion \[real, dual\], nx8
%[text] invDq: inverse dual quaternion \[real, dual\]，nx8
%[text] ## note
%[text] 
%[text] ## references 
%[text] Sveier, A., & Egeland, O. (2020). Dual Quaternion Particle Filtering for Pose Estimation. IEEE Transactions on Control Systems Technology, 1-14.
%[text] ## revisions
%[text] 20220526  y.yoshimura
%[text] See also qMult, dqConj.
function invDq = dqInv(scalar, dq)
% arguments
%     scalar (1,1) {mustBeMember(scalar, [0, 4])}
%     dq (:,8) {mustBeNumeric}
% end

dqr = dq(:,1:4); % real part
dqd = dq(:,5:8); % dual part

invDqr = qInv(scalar, dqr);
invDqd = -qMult(scalar, 0, qMult(scalar, 0, invDqr, dqd), invDqr);

invDq = [invDqr invDqd];


end

%[appendix]{"version":"1.0"}
%---
