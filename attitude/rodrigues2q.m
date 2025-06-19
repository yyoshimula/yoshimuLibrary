%[text] # calculating quaternions from Rodrigues parameters
%[text] `scalar:` specify the definition of the quaternion 
%[text] `scalar == 0`
%[text] ${\\bf q} = \[q\_0,q\_1,q\_2,q\_3\]^T=\[\\cos(\\theta/2), {\\bf e}^T\\sin(\\theta/2)\]^T$
%[text] `scalar == 4`
%[text] ${\\bf q} = \[q\_1,q\_2,q\_3, q\_4\]^T=\[{\\bf e}^T\\sin(\\theta/2), \\cos(\\theta/2)\]^T$
%[text] `q`: quaternions, nx4 matrix
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20150101  y.yoshimura, y.yoshimula@gmail.com
%[text] See also q2zyx.
function q = rodrigues2q(scalar, rod)
% arguments
%     scalar (1,1) {mustBeMember(scalar, [0, 4])}
%     rod (:,3) {mustBeNumeric}
% end

f = 1;
a = 0;

q = grp2q(scalar, f, a, rod);

end


%[appendix]{"version":"1.0"}
%---
