%[text] # calculating quaternions from direction cosine matrix (DCM)
%[text] 回転行列からquaternionを計算
%[text] ただしこの方法だとquaternionが不連続になるので注意．連続性も確保したい場合は dcm2qC.mlxを使う．
%[text] ## inputs
%[text] `scalar:` specify the definition of the quaternion 
%[text] `scalar == 0`
%[text] ${\\bf q} = \[q\_0,q\_1,q\_2,q\_3\]^T=\[\\cos(\\theta/2), {\\bf e}^T\\sin(\\theta/2)\]^T$
%[text] `scalar == 4`
%[text] ${\\bf q} = \[q\_1,q\_2,q\_3, q\_4\]^T=\[{\\bf e}^T\\sin(\\theta/2), \\cos(\\theta/2)\]^T$\]
%[text] R, rotation matrix, 3x3 matrix 
%[text] ## output
%[text] `q`: quaternions, 1x4 vector 
%[text] ## note
%[text] NA
%[text] ## references 
%[text]  Markley, F. L., “Unit Quaternion from Rotation Matrix,” Journal of Guidance Control, and Dynamics, vol. 31, Mar. 2008, pp. 440-442.
%[text] ## revisions
%[text] 20150101  y.yoshimura, y.yoshimula@gmail.com
%[text] See also q2DCM.
function q = dcm2q(scalar, R)
% arguments
%     scalar (1,1) {mustBeMember(scalar, [0, 4])}
%     R (3,3) {mustBeNumeric}
% end

trR =  R(1,1) + R(2,2) + R(3,3); % trace
tmp = [R(1,1)
    R(2,2)
    R(3,3)
    trR];

[~, iq] = max(tmp);
switch iq

    case 1
        x1 = [1 + 2*R(1,1) - trR
            R(1,2) + R(2,1)
            R(1,3) + R(3,1)
            R(2,3) - R(3,2)];

        q = x1 ./ norm(x1);

    case 2
        x2 = [R(2,1) + R(1,2)
            1 + 2 * R(2,2) - trR
            R(2,3) + R(3,2)
            R(3,1) - R(1,3)];

        q = x2 ./ norm(x2);
    case 3
        x3 = [R(3,1) + R(1,3)
            R(3,2) + R(2,3)
            1 + 2 * R(3,3) - trR
            R(1,2) - R(2,1)];

        q = x3 ./ norm(x3);
    case 4
        x4 = [R(2,3) - R(3,2)
            R(3,1) - R(1,3)
            R(1,2) - R(2,1)
            1 + trR];

        q = x4 ./ norm(x4);
end

q = q';

if scalar == 0
    q = [q(4), q(1), q(2), q(3)];
end

end

%[appendix]{"version":"1.0"}
%---
