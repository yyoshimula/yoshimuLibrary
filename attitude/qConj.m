%[text] # quaternion conjugate
%[text] `scalar:` specify the definition of the quaternion 
%[text] `scalar == 0`
%[text] ${\\bf q} = \[q\_0,q\_1,q\_2,q\_3\]^T=\[\\cos(\\theta/2), {\\bf e}^T\\sin(\\theta/2)\]^T$
%[text] `scalar == 4`
%[text] ${\\bf q} = \[q\_1,q\_2,q\_3, q\_4\]^T=\[{\\bf e}^T\\sin(\\theta/2), \\cos(\\theta/2)\]^T$
%[text] `q`: quaternions, 1x4 vector 
%[text] `qC`, quaternion conjugate, nx4 matrix
%[text] ${\\bf q}^\\dagger=\[-{\\bf q}\_v^T, q\_4\]^T$ 
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20150101  y.yoshimura, y.yoshimula@gmail.com
%[text] See also qInv.
function qC = qConj(scalar, q)
% arguments
%     scalar (1,1) {mustBeMember(scalar, [0, 4])}
%     q (:,4)
% end

if scalar == 0
    qC = [q(:,1), -q(:,2:4)];
    
elseif scalar == 4
    qC = [-q(:,1:3) q(:,4)];
    
else
    error('definition of quaternions is unclear');
    
end

end

%[appendix]{"version":"1.0"}
%---
