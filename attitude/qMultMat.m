%[text] # Quaternion multiplication matrix
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
%[text] `q`: quaternions, 1x4 matrix
%[text] qMat: 4x4 matrix 
%[text] ## note
%[text] NA
%[text] ## references 
%[text] Markley, F. L., and Crassidis, J. L., Fundamentals of Spacecraft Attitude Determination and Control, New York, NY: Springer, 2014. 
%[text] ## revisions
%[text] 20231219  y.yoshimura, y.yoshimula@gmail.com
%[text] See also qCon, qInv.
function qMat = qMultMat(scalar, def, q)
% arguments
%     scalar (1,1) {mustBeMember(scalar, [0, 4])}
%     def (1,1) {mustBeMember(def, [0, 1])}
%     q (:,4) {mustBeNumeric}
% end

if (scalar == 0)

    q0 = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);
    qv = [q1; q2; q3];

    if def == 0
        qMat = [q0, -qv'
            qv, q0 .* eye(3) + skew(qv)];
    elseif def == 1
        qMat=  [q0, -qv'
            qv, q0 .* eye(3) - skew(qv)];
    else
        error('definition of quaternion multiplication is unclear');
    end


elseif scalar == 4
    q4 = q(4);
    q1 = q(1);
    q2 = q(2);
    q3 = q(3);
    qv = [q1; q2; q3];

    if def == 0
        qMat = [q4.*eye(3)+skew(qv), qv
        -qv' q4];
    elseif def == 1
        qMat = [q4.*eye(3)-skew(qv), qv
        -qv' q4];
    else
        error('definition of quaternion multiplication is unclear');
    end

else
    error('definition of quaternions is unclear');

end

%[appendix]{"version":"1.0"}
%---
