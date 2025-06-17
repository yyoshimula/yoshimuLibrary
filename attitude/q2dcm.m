%[text] # calculating directional cosine matrix (DCM) from quaternions 
%[text] quaternionから回転行列を計算
%[text] `scalar:` specify the definition of the quaternion 
%[text] `scalar == 0`
%[text] ${\\bf q} = \[q\_0,q\_1,q\_2,q\_3\]^T=\[\\cos(\\theta/2), {\\bf e}^T\\sin(\\theta/2)\]^T$
%[text] `scalar == 4`
%[text] ${\\bf q} = \[q\_1,q\_2,q\_3, q\_4\]^T=\[{\\bf e}^T\\sin(\\theta/2), \\cos(\\theta/2)\]^T$
%[text] `q`: quaternions, 1x4 vector 
%[text] `DCM`, rotation matrix, 3x3 matrix 
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20150101  y.yoshimura, y.yoshimula@gmail.com
%[text] See also q2zyx.
function DCM = q2dcm(scalar, q)
% arguments    
%     scalar (1,1) {mustBeMember(scalar, [0, 4])}
%     q (:,4) {mustBeNumeric}
% end

qNorm = sqrt(q(1)^2 + q(2)^2 + q(3)^2 + q(4)^2);
qTmp = q ./ qNorm; % normalize

% calculate using q4 definition
if (scalar == 0)
    q4 = qTmp(1);
    q1 = qTmp(2);
    q2 = qTmp(3);
    q3 = qTmp(4);
        
elseif (scalar == 4)
    q1 = qTmp(1);
    q2 = qTmp(2);
    q3 = qTmp(3);
    q4 = qTmp(4);
    
end

DCM = [q1^2 - q2^2 - q3^2 + q4^2, 2.0 * (q1 * q2 + q4 * q3), 2.0 * (q1 * q3 - q4 * q2)
        2.0 * (q1 * q2 - q4 * q3), - q1^2 + q2^2 - q3^2 + q4^2, 2.0 * (q2 * q3 + q4 * q1)
        2.0 * (q1 * q3 + q4 * q2), 2.0 * (q2 * q3 - q4 * q1), - q1^2 - q2^2 + q3^2 + q4^2 ];
    

end

%[appendix]{"version":"1.0"}
%---
