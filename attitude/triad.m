%[text] # calculate diretional cosine matrix (DCM) using TRIAD method
%[text] ## input
%[text] `v1_i, v2_i`: reference vectors expressed with inertial frame, nx3 matrix
%[text] `w1_b, w2_b`: reference vectors expressed with body-fixed frame, nx3 matrix
%[text] ## output
%[text] R: directional cosine matrix (rotation matrix) from inertial frame to body-fixed frame, 3 x 3 x n matrix
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20220201  y.yoshimura, y.yoshimula@gmail.com
%[text] See also q2dcm, dcm2q.
function R = triad(v1_i, v2_i, w1_b, w2_b)
% arguments
%     v1_i (:,3) {mustBeNumeric}
%     v2_i (:,3) {mustBeNumeric}
%     w1_b (:,3) {mustBeNumeric}
%     w2_b (:,3) {mustBeNumeric}
% end

% normalization
v1_i = v1_i ./ vecnorm(v1_i, 2, 2);
v2_i = v2_i ./ vecnorm(v2_i, 2, 2);
w1_b = w1_b ./ vecnorm(w1_b, 2, 2);
w2_b = w2_b ./ vecnorm(w2_b, 2, 2);

n = size(v1_i, 1);

% pre-allocation
R = zeros(3,3,n);

for i = 1:n
    r1 = v1_i(i,:);
    r2 = cross(v1_i(i,:), v2_i(i,:));
    r3 = cross(r1, r2);

    s1 = w1_b(i,:);
    s2 = cross(w1_b(i,:), w2_b(i,:));
    s3 = cross(s1, s2);

    R(:,:,i) = [s1' s2' s3'] * [r1' r2' r3']';
end


end

%[appendix]{"version":"1.0"}
%---
