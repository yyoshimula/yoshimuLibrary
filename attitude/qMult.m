%[text] # Quaternion multiplication
%[text] `scalar:` specifies the definition of the quaternion 
%[text] `scalar == 0`
%[text] ${\\bf q} = \[q\_0,q\_1,q\_2,q\_3\]^T=\[\\cos(\\theta/2), {\\bf e}^T\\sin(\\theta/2)\]^T$
%[text] `scalar == 4`
%[text] ${\\bf q} = \[q\_1,q\_2,q\_3, q\_4\]^T=\[{\\bf e}^T\\sin(\\theta/2), \\cos(\\theta/2)\]^T$
%[text] `def`: specifies the definition of the quaternion multiplication
%[text] `def == 0`
%[text] ${\\bf q}\\odot{\\bf p}= \\left\[\\begin{array}{cc}\n     q\_{4}\\bar{{\\bf{p}}}+p\_{4}\\bar{{\\bf{q}}} + \\bar{{\\bf{q}}}\\times\\bar{{\\bf{p}}}\\\\\nq\_{4}p\_{4}-\\bar{{\\bf{q}}}^{T}\\bar{{\\bf{p}}}\n\\end{array}\\right\]$
%[text] `def == 1`
%[text] ${\\bf q}\\otimes{\\bf p}= \\left\[\\begin{array}{cc}\n     q\_{4}\\bar{{\\bf{p}}}+p\_{4}\\bar{{\\bf{q}}} - \\bar{{\\bf{q}}}\\times\\bar{{\\bf{p}}}\\\\\nq\_{4}p\_{4}-\\bar{{\\bf{q}}}^{T}\\bar{{\\bf{p}}}\n\\end{array}\\right\]$
%[text] `q`: quaternions, nx4 matrix
%[text] `p`: quaternions, nx4 matrix 
%[text] `output`: quaternions, nx4 matrix 
%[text] ## note
%[text] NA
%[text] ## references 
%[text] Markley, F. L., and Crassidis, J. L., Fundamentals of Spacecraft Attitude Determination and Control, New York, NY: Springer, 2014. 
%[text] ## revisions
%[text] 20150101  y.yoshimura, y.yoshimula@gmail.com
%[text] See also qCon, qInv.
function output = qMult(scalar, def, q, p)
% arguments
%     scalar (1,1) {mustBeMember(scalar, [0, 4])}
%     def (1,1) {mustBeMember(def, [0, 1])}
%     q (:,4)
%     p (:,4)
% end


% 入力のサイズを取得
nQ = size(q, 1);
output = zeros(nQ, 4);

if (scalar == 0)

    q0 = q(:,1);
    q1 = q(:,2);
    q2 = q(:,3);
    q3 = q(:,4);

    p0 = p(:,1);
    p1 = p(:,2);
    p2 = p(:,3);
    p3 = p(:,4);

    % スカラー部の計算
    scalarPart = q0 .* p0 - (q1 .* p1 + q2 .* p2 + q3 .* p3);

    % ベクトル部の計算（cross積を直接展開）
    % cross(qv, pv)
    cross1 = q2 .* p3 - q3 .* p2;
    cross2 = q3 .* p1 - q1 .* p3;
    cross3 = q1 .* p2 - q2 .* p1;

    if def == 0

        output = [scalarPart, ...
            q0 .* p1 + p0 .* q1 + cross1, ...
            q0 .* p2 + p0 .* q2 + cross2, ...
            q0 .* p3 + p0 .* q3 + cross3];

    else % def == 1
        output = [scalarPart, ...
            q0 .* p1 + p0 .* q1 - cross1, ...
            q0 .* p2 + p0 .* q2 - cross2, ...
            q0 .* p3 + p0 .* q3 - cross3];
    end

else % (scalar == 4)
    q1 = q(:,1);
    q2 = q(:,2);
    q3 = q(:,3);
    q4 = q(:,4);

    p1 = p(:,1);
    p2 = p(:,2);
    p3 = p(:,3);
    p4 = p(:,4);

    scalarPart = q4 .* p4 - (q1 .* p1 + q2 .* p2 + q3 .* p3);

    % （cross積を直接展開 = cross(qv, pv) ）
    cross1 = p3 .* q2 - p2 .* q3;
    cross2 = p1 .* q3 - p3 .* q1;
    cross3 = p2 .* q1 - p1 .* q2;
    if def == 0
        output = [q4 .* p1 + p4 .* q1 + cross1, ...
            q4 .* p2 + p4 .* q2 + cross2, ...
            q4 .* p3 + p4 .* q3 + cross3, ...
            scalarPart];
    else % def == 1
        output = [q4 .* p1 + p4 .* q1 - cross1, ...
            q4 .* p2 + p4 .* q2 - cross2, ...
            q4 .* p3 + p4 .* q3 - cross3, ...
            scalarPart];
    end

end
end

%[appendix]{"version":"1.0"}
%---
