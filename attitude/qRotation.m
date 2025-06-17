%[text] # rotating coordiante frames with quaternions
%[text] ベクトルを回転後の座標系で表示（ベクトルを回転させているわけではない）
%[text] `scalar:` specify the definition of the quaternion 
%[text] `scalar == 0`
%[text] ${\\bf q} = \[q\_0,q\_1,q\_2,q\_3\]^T=\[\\cos(\\theta/2), {\\bf e}^T\\sin(\\theta/2)\]^T$
%[text] `scalar == 4`
%[text] ${\\bf q} = \[q\_1,q\_2,q\_3, q\_4\]^T=\[{\\bf e}^T\\sin(\\theta/2), \\cos(\\theta/2)\]^T$
%[text] `r:` vector, nx3 vector
%[text] `q`: quaternions, nx4 vector 
%[text] `rb`: vector expressed with the rotated coordinate frame, 
%[text] ## note
%[text] ${\\bf{r}}\_{b} = {\\bf{q}} \\otimes {\\bf r} \\otimes {\\bf q}^{-1} = {\\bf q}^{-1}\\odot {\\bf r}\\odot {\\bf q}$
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20150101  y.yoshimura, y.yoshimula@gmail.com
%[text] See also qInv.
function rb = qRotation(scalar, r, q)
% arguments    
%     scalar (1,1) {mustBeMember(scalar, [0, 4])}
%     r (:,3) {mustBeNumeric}
%     q (:,4) {mustBeNumeric}
% end

n = size(q,1);
rTmp = [zeros(n,1) r];

% q0がスカラーの定義かつ \odot のクォータニオン積の定義で計算する
q = (scalar == 0) .* q ...
    + (scalar == 4) .* [q(:,4), q(:,1), q(:,2), q(:,3)];

qInverse = qInv(0, q);

tmp  = qMult(0, 0, qInverse, rTmp);
tmp2 = qMult(0, 0, tmp, q);

rb = tmp2(:,2:4);

end

%[appendix]{"version":"1.0"}
%---
