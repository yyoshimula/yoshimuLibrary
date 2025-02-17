% quaternion inverse
% scalar: specify the definition of the quaternion
% scalar == 0
%
% scalar == 4
%
% q: quaternions, 1x4 vector
% qInv, quaternion inverse, nx4 matrix
%     where  is quaternion conjugate
% note
% NA
% references
% NA
% revisions
% 20150101  y.yoshimura, y.yoshimula@gmail.com
% See also qConj.

function qInv = qInv_(scalar, q)

if scalar == 0
    qInv = [q(:,1) -q(:,2) -q(:,3) -q(:,4)] ./ vecnorm(q,2,2).^2;

elseif scalar == 4
    qInv = [-q(:,1) -q(:,2) -q(:,3) q(:,4)] ./ vecnorm(q,2,2).^2;
else
    error('invalid qwuaternion definition')
end

end
