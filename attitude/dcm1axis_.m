% direction cosine matrix (DCM) about a single axis
% 1軸周りの回転行列を計算
% axis: rotation axis, 1 == x-axis, 2 == y-axis, 3 == z-axis
% phi: rotation angle, rad, scalar
% R, rotation matrix, 3x3 matrix 
% note
% NA
% references 
% NA
% revisions
% 20150101  y.yoshimura, y.yoshimula@gmail.com
% See also q2DCM, dcm2q.

function R = dcm1axis_(axis, phi)

cp = cos(phi);
sp = sin(phi);

R = (axis == 1) .* [1 0 0
    0 cp sp
    0 -sp cp]...
    + (axis == 2) .* [cp 0 -sp
    0 1 0
    sp 0 cp]...
    + (axis == 3) .* [cp sp 0
    -sp cp 0
    0 0 1];

end
