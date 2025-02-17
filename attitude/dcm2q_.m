function q = dcm2q_(scalar, R)
% calculating quaternions from direction cosine matrix (DCM)
% 回転行列からquaternionを計算
% ただしこの方法だとquaternionが不連続になるので注意．連続性も確保したい場合は dcm2qC.mlxを使う．
% inputs
% scalar: specify the definition of the quaternion
% scalar == 0
%
% scalar == 4
% ]
% R, rotation matrix, 3x3 matrix
% output
% q: quaternions, 1x4 vector
% note
% NA
% references
%  Markley, F. L., “Unit Quaternion from Rotation Matrix,” Journal of Guidance Control, and Dynamics, vol. 31, Mar. 2008, pp. 440-442.
% revisions
% 20150101  y.yoshimura, y.yoshimula@gmail.com
% See also q2DCM.

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

q = (scalar == 0 ) .* [q(4) q(1) q(2) q(3)]...
    + (scalar == 4) .* q';

end
