%[text] # Cayley transformation
%[text] ケイリー変換
%[text] ${\\rm cay}(\\bf{u}) = (1+\\bf{u})\\odot (1-\\bf{u})^{-1}$ 
%[text] scalar: 0: quaternionのスカラ部q0としてq=\[q0 q1 q2 q3\]という定義
%[text]             4: quaternionのスカラ部q4としてq=\[q1 q2 q3 q4\]という定義
%[text] u: 任意のベクトル，nx3
%[text] q: quaternion, nx4
%[text] ## note
%[text] 定義式における1-uをまとめてquaternion．つまり \[$\[{\\bf u}^T, 1\]^T$のこと．
%[text] ## references 
%[text] Sveier, A., & Egeland, O. (2020). Dual Quaternion Particle Filtering for Pose Estimation. IEEE Transactions on Control Systems Technology, 1-14.
%[text] ## revisions
%[text] 20220526  y.yoshimura
%[text] See also qMult
function q = cay(scalar, u)
% arguments
%     scalar (1,1) {mustBeMember(scalar, [0 4])}
%     u (:,3) {mustBeNumeric}
% end

n = size(u,1);

% u(4) = scalar partとして計算
if scalar == 0
    uQ1 = [ones(n,1), u];
elseif scalar == 4
    uQ1 = [u, ones(n,1)];
else
    error('quaternion definition is unclear')
end


uQ2 = qConj(4, uQ1);

% u(4) = scalar partとして計算
q = qMult(4, 0, uQ1, qInv(scalar, uQ2));

if scalar == 0
    q = [q(:,4), q(:,1:3)];
elseif scalar == 4
    % as is
else
    error('quaternion definition is unclear')
end

end


%[appendix]{"version":"1.0"}
%---
