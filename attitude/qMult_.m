% Quaternion multiplication
% scalar: specifies the definition of the quaternion
% scalar == 0
%
% scalar == 4
%
% def: specifies the definition of the quaternion multiplication
% def == 0
%
% def == 1
%
% q: quaternions, nx4 matrix
% p: quaternions, nx4 matrix
% output: quaternions, nx4 matrix
% note
% NA
% references
% Markley, F. L., and Crassidis, J. L., Fundamentals of Spacecraft Attitude Determination and Control, New York, NY: Springer, 2014.
% revisions
% 20150101  y.yoshimura, y.yoshimula@gmail.com
% See also qCon, qInv.

function output = qMult_(scalar, def, q, p)

if (scalar == 0)

    q0 = q(:,1);
    qv = q(:,2:4);
    q1 = qv(:,1);
    q2 = qv(:,2);
    q3 = qv(:,3);

    p0 = p(:,1);
    pv = p(:,2:4);
    p1 = pv(:,1);
    p2 = pv(:,2);
    p3 = pv(:,3);

    if def == 0
        output = [q0 .* p0 - (q1 .* p1 + q2 .* p2 + q3 .* p3), ...
            q0 .* pv + p0 .* qv + cross(qv, pv)];

    else
        output = [q0 .* p0  - (q1 .* p1 + q2 .* p2 + q3 .* p3), ...
            q0 .* pv + p0 .* qv - cross(qv, pv)];
    end

else % (scalar == 4)
    qv = q(:,1:3);
    q1 = qv(:,1);
    q2 = qv(:,2);
    q3 = qv(:,3);
    q4 = q(:,4);

    pv = p(:,1:3);
    p1 = pv(:,1);
    p2 = pv(:,2);
    p3 = pv(:,3);
    p4 = p(:,4);

    if def == 0
        output = [q4 .* pv + p4 .* qv + cross(qv,pv), ...
            q4 .* p4 - (q1 .* p1 + q2 .* p2 + q3 .* p3)];
    else
        output = [q4 .* pv + p4 .* qv - cross(qv,pv), ...
            q4 .* p4 - (q1 .* p1 + q2 .* p2 + q3 .* p3)];
    end

end