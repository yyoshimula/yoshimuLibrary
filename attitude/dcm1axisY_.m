% direction cosine matrix (DCM) about a single axis
% 1軸周りの回転行列を計算
% phi: rotation angle, rad, scalar
% R, rotation matrix, 3x3 matrix 
% note
% NA
% references 
% NA
% revisions
% 20150101  y.yoshimura, y.yoshimula@gmail.com
% See also q2DCM, dcm2q.

function R = dcm1axisY_(phi)

cp = cos(phi);
sp = sin(phi);

R = [cp 0 -sp
    0 1 0
    sp 0 cp];

end
