%[text] # Inverse Cayley transformation
%[text] 逆ケイリー変換
%[text] ${\\rm cay}^{-1}(\\bf{q}) = (\\bf{q}-1)\\odot (\\bf{q}+1)^{-1}$ 
%[text] scalar: 0: quaternionのスカラ部q0としてq=\[q0 q1 q2 q3\]という定義
%[text]             4: quaternionのスカラ部q4としてq=\[q1 q2 q3 q4\]という定義
%[text] q: quaternion, nx4
%[text] u: 任意のベクトル，nx3
%[text] ## note
%[text] 定義式における1はunit quaternion．つまり \[0 0 0 1\]のこと．
%[text] ## references 
%[text] Sveier, A., & Egeland, O. (2020). Dual Quaternion Particle Filtering for Pose Estimation. IEEE Transactions on Control Systems Technology, 1-14.
%[text] ## revisions
%[text] 20220526  y.yoshimura
%[text] See also qMult
function u = cayInv(scalar, q)
% arguments
%     scalar (1,1) {mustBeMember(scalar, [0 4])}
%     q (:,4) {mustBeNumeric}
% end

n = size(q,1);
% q4 = scalar partとして計算
if scalar == 0
    q = [q(:,2:4), q(:,1)];
elseif scalar == 4
    % as is
else
    error('quaternion definition is unclear')
end

qUnit = [zeros(n,3), ones(n,1)];

tmp = qMult(4, 0, q-qUnit, qInv(4,q+qUnit));

u = tmp(:,1:3);

end


%[appendix]{"version":"1.0"}
%---
