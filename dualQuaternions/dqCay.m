%[text] # dual Cayley transformation
%[text] dual vectorのケイリー変換
%[text] ${\\rm cay}(\\bf{\\tilde{u}}) = (1+\\bf{\\tilde{u}})\\odot (1-\\bf{\\tilde{u}})^{-1}$ 
%[text] ↑の1+はreal partのスカラ部を1にするという意味．= $1+\\epsilon 0$を足すということ
%[text] scalar: 0: quaternionのスカラ部q0としてq=\[q0 q1 q2 q3\]という定義
%[text]             4: quaternionのスカラ部q4としてq=\[q1 q2 q3 q4\]という定義
%[text] `dv`: 任意のdualベクトル（なのでスカラ部は0），nx8
%[text] `dq`: dual quaternion, nx8
%[text] ## note
%[text] ## references 
%[text] Sveier, A., & Egeland, O. (2020). Dual Quaternion Particle Filtering for Pose Estimation. IEEE Transactions on Control Systems Technology, 1-14.
%[text] ## revisions
%[text] 20240111  y.yoshimura
%[text] See also cay.
function dq = dqCay(scalar, dv)
% arguments
%     scalar (1,1) {mustBeMember(scalar, [0 4])}
%     dv (:,8) {mustBeNumeric}
% end

n = size(dv,1);

% real partはスカラ部を1に
if scalar == 0
    dvReal = [ones(n,1), dv(:,2:4)];
    tmp = dv(:,2:4); % for real part calculation
else
    dvReal = [dv(:,1:3), ones(n,1)];
    tmp = dv(:,1:3);
end

% dual part 一応、スカラ部をゼロに．
% dvDual = (scalar == 0) .* [zeros(n,1), dv(:,6:8)] ...
%     + (scalar == 4) .* [dv(:,5:7), zeros(n,1)];

dvDual = dv(:,5:8);
%[text] ## Cayley transformation
% real part
dqReal = cay(scalar, tmp);

% dual part
dqDual = 2 * qMult(scalar, 0, qMult(scalar, 0, dvReal, dvDual), dvReal);
dqDual = dqDual ./ (vecnorm(dvReal,2,2).^2).^2;

dq = [dqReal, dqDual];

end


%[appendix]{"version":"1.0"}
%---
