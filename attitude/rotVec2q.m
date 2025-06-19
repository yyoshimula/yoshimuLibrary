%[text] # calculating quaternionf from rotation vectors
%[text] `scalar:` specify the definition of the quaternion 
%[text] `scalar == 0`
%[text] ${\\bf q} = \[q\_0,q\_1,q\_2,q\_3\]^T=\[\\cos(\\theta/2), {\\bf e}^T\\sin(\\theta/2)\]^T$
%[text] `scalar == 4`
%[text] ${\\bf q} = \[q\_1,q\_2,q\_3, q\_4\]^T=\[{\\bf e}^T\\sin(\\theta/2), \\cos(\\theta/2)\]^T$
%[text] `q`: quaternions, nx4 matrix
%[text] `rv:` $\\theta \\bf e$
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20230614  y.yoshimura, y.yoshimula@gmail.com
%[text] See also q2zyx.
function q = rotVec2q(scalar, rv)
% arguments
%     scalar (1,1) {mustBeMember(scalar, [0, 4])}
%     rv (:,3) {mustBeNumeric}
% end

theta = vecnorm(rv, 2, 2);
eAxis = rv ./ theta;

if scalar == 0
    q = [cos(theta./2), eAxis.*sin(theta./2)];
elseif scalar == 4
    q = [eAxis.*sin(theta./2), cos(theta./2)];
else
    error('quaternion definition is unlcear')
end
    

end



%[appendix]{"version":"1.0"}
%---
