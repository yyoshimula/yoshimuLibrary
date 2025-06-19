%[text] # quaternion around a single axis
%[text] 1軸周りのquaternionを計算
%[text] `scalar:` specify the definition of the quaternion 
%[text] `scalar == 0`
%[text] ${\\bf q} = \[q\_0,q\_1,q\_2,q\_3\]^T=\[\\cos(\\theta/2), {\\bf e}^T\\sin(\\theta/2)\]^T$
%[text] `scalar == 4`
%[text] ${\\bf q} = \[q\_1,q\_2,q\_3, q\_4\]^T=\[{\\bf e}^T\\sin(\\theta/2), \\cos(\\theta/2)\]^T$
%[text] `axis`: rotation axis, unit vector, nx3 matrix
%[text] `theta`: rotation angle, rad, scalar, nx1 vector
%[text] q, quaternion, nx4 matrix 
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20150101  y.yoshimura, y.yoshimula@gmail.com
%[text] See also dcm1axis.
function q = q1axis(scalar, axis, theta)
% arguments
%     scalar (1,1) {mustBeMember(scalar,[0, 4])}
%     axis (:,3) {mustBeNumeric}
%     theta (:,1) {mustBeNumeric}
% end

% normalize
axis = axis ./ vecnorm(axis,2,2);

%[text] ${\\bf q} = \[q\_1,q\_2,q\_3, q\_4\]^T=\[{\\bf e}^T\\sin(\\theta/2), \\cos(\\theta/2)\]^T$で一旦計算
q = [axis .* sin(theta./2), cos(theta./2)];

if scalar == 0
    q = [q(:,4) q(:,1) q(:,2) q(:,3)];
end

end

%[appendix]{"version":"1.0"}
%---
