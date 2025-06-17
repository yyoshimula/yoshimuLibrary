%[text] # quaternion inverse
%[text] `scalar:` specify the definition of the quaternion 
%[text] `scalar == 0`
%[text] ${\\bf q} = \[q\_0,q\_1,q\_2,q\_3\]^T=\[\\cos(\\theta/2), {\\bf e}^T\\sin(\\theta/2)\]^T$
%[text] `scalar == 4`
%[text] ${\\bf q} = \[q\_1,q\_2,q\_3, q\_4\]^T=\[{\\bf e}^T\\sin(\\theta/2), \\cos(\\theta/2)\]^T$
%[text] `q`: quaternions, 1x4 vector 
%[text] `qInv`, quaternion inverse, nx4 matrix
%[text] ${\\bf q}^{-1} =\\frac{{\\bf q}^\\ast}{\\|{\\bf q}\\|^2}$    where ${\\bf q}^\\ast\n$ is quaternion conjugate
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20150101  y.yoshimura, y.yoshimula@gmail.com
%[text] See also qConj.
function qInv = qInv(scalar, q)
% arguments    
%     scalar (1,1) {mustBeMember(scalar, [0, 4])}
%     q (:,4) {mustBeNumeric}
% end

qInv = (scalar == 0 ) .* [q(:,1) -q(:,2) -q(:,3) -q(:,4)] ./ vecnorm(q,2,2).^2 ...
    + (scalar == 4) .* [-q(:,1) -q(:,2) -q(:,3) q(:,4)] ./ vecnorm(q,2,2).^2;

end

%[appendix]{"version":"1.0"}
%---
