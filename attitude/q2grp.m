%[text] # converting quaternions to generalized Rodrigues parameters
%[text] `scalar == 0`
%[text] ${\\bf q} = \[q\_0,q\_1,q\_2,q\_3\]^T=\[\\cos(\\theta/2), {\\bf e}^T\\sin(\\theta/2)\]^T$
%[text] `scalar == 4`
%[text] ${\\bf q} = \[q\_1,q\_2,q\_3, q\_4\]^T=\[{\\bf e}^T\\sin(\\theta/2), \\cos(\\theta/2)\]^T$
%[text] `f:` scaling parameter for Rodrigues parameters, scalar
%[text] `a:` parameter for Rodrigues parameters, scalar 
%[text] `q`, quaternions, nx4 matrix
%[text] `p`, generalized Rodrigues parameters, nx3 matrix
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20190121  y.yoshimura, y.yoshimula@gmail.com
%[text] See also GRP2q.
function p = q2grp(scalar, f, a, q)
% arguments
%     scalar (1,1) {mustBeMember(scalar, [0, 4])}
%     f (1,1) {mustBeNumeric}
%     a (1,1) {mustBeNumeric}
%     q (:,4) {mustBeNumeric}
% end

% q4 = scalar partとして計算．
if scalar == 0
q = [q(:,2:4), q(:,1)];
end

p = f .* sign(q(:,4)) .* q(:,1:3) ./ (a + abs(q(:,4)));

end

%[appendix]{"version":"1.0"}
%---
