%[text] # spherical linear interpolation (slerp)
%[text] quaternionの線形補間
%[text] t: normalized time, i.e., 0 \<= t \<= 1 
%[text] `scalar:` specify the definition of the quaternion 
%[text] `scalar == 0`
%[text] ${\\bf q} = \[q\_0,q\_1,q\_2,q\_3\]^T=\[\\cos(\\theta/2), {\\bf e}^T\\sin(\\theta/2)\]^T$
%[text] `scalar == 4`
%[text] ${\\bf q} = \[q\_1,q\_2,q\_3, q\_4\]^T=\[{\\bf e}^T\\sin(\\theta/2), \\cos(\\theta/2)\]^T$
%[text] `q`: quaternions, 1x4 vector 
%[text] `qt`: interpolated quaternions, 1x4 vector 
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20150101  y.yoshimura, y.yoshimula@gmail.com
%[text] See also sclerp.
function qt = slerp(t, scalar, q1, q2)

t = t(:); % column vector

qtmp = qMult(scalar, 0 , qInv(scalar,q1), q2);

if scalar == 0
    eTheta = acos(qtmp(1));
    eTheta = eTheta * 2;
    eAxis = qtmp(2:4);

elseif scalar == 4
    eTheta = acos(qtmp(4)); % Euler angle of rotation
    eTheta = eTheta * 2;
    eAxis = qtmp(1:3);

else
    error('quaternion definition is unclear')
end

eAxis = eAxis ./ sin(eTheta/2); % Euler axis of rotation

q1tmp = kron(ones(length(t),1), q1); % Nx4

if scalar == 0
    qt = qMult(scalar, 0, q1tmp, [cos(t.*eTheta/2) eAxis.*sin(t.*eTheta/2)]);

elseif scalar == 4
    qt = qMult(scalar, 0, q1tmp, [eAxis.*sin(t.*eTheta/2) cos(t.*eTheta/2)]);
    
else
    error('quaternion definition is unclear')

end

end

%[appendix]{"version":"1.0"}
%---
