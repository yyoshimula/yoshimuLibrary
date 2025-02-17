% rotating coordiante frames with quaternions
% ベクトルを回転後の座標系で表示（ベクトルを回転させているわけではない）
% scalar: specify the definition of the quaternion 
% scalar == 0
% 
% scalar == 4
% 
% r: vector, nx3 vector
% q: quaternions, nx4 vector 
% rb: vector expressed with the rotated coordinate frame, 
% note
% 
% references 
% NA
% revisions
% 20150101  y.yoshimura, y.yoshimula@gmail.com
% See also qInv.

function rb = qRotation_(scalar, r, q)

n = max(size(q,1), size(r,1));
rTmp = [zeros(n,1) r];

% q0がスカラーの定義かつ \odot のクォータニオン積の定義で計算する
q = (scalar == 0) .* q ...
    + (scalar == 4) .* [q(:,4), q(:,1), q(:,2), q(:,3)];

qInverse = qInv_(0, q);

tmp  = qMult_(0, 0, qInverse, rTmp);
tmp2 = qMult_(0, 0, tmp, q);

rb = tmp2(:,2:4);

end
