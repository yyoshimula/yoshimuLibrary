%[text] # direction cosine matrix (DCM) about Z axis
%[text] x軸周りの回転行列を計算
%[text] `phi`: rotation angle, rad, scalar
%[text] R, rotation matrix, 3x3 matrix 
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20241216  y.yoshimura, y.yoshimula@gmail.com
%[text] See also q2DCM, dcm1axis.
function R = dcm1axisZ(phi)
% arguments        
%     phi (1,1) 
% end

cp = cos(phi);
sp = sin(phi);

R = [cp sp 0
    -sp cp 0
    0 0 1];

end

%[appendix]{"version":"1.0"}
%---
