%[text] # direction cosine matrix (DCM) about a single axis
%[text] 1軸周りの回転行列を計算
%[text] `axis`: rotation axis, 1 == x-axis, 2 == y-axis, 3 == z-axis
%[text] `phi`: rotation angle, rad, scalar
%[text] R, rotation matrix, 3x3 matrix 
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20150101  y.yoshimura, y.yoshimula@gmail.com
%[text] See also q2DCM, dcm2q.
function R = dcm1axis(axis, phi)
% arguments
%     axis (1,1) {mustBeMember(axis, [1, 2, 3])}
%     phi (1,1)
% end

cp = cos(phi);
sp = sin(phi);

if axis == 1
    R = [1 0 0
        0 cp sp
        0 -sp cp];
elseif axis == 2
    R = [cp 0 -sp
        0 1 0
        sp 0 cp];
else
    R = [cp sp 0
        -sp cp 0
        0 0 1];
end
end

%[appendix]{"version":"1.0"}
%---
