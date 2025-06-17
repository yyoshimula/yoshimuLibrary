%[text] # dual quaternion multiplication
%[text] `scalar:` specify the definition of the quaternion 
%[text] `scalar == 0`
%[text] ${\\bf q} = \[q\_0,q\_1,q\_2,q\_3\]^T=\[\\cos(\\theta/2), {\\bf e}^T\\sin(\\theta/2)\]^T$
%[text] `scalar == 4`
%[text] ${\\bf q} = \[q\_1,q\_2,q\_3, q\_4\]^T=\[{\\bf e}^T\\sin(\\theta/2), \\cos(\\theta/2)\]^T$
%[text] `def`: specifies the definition of the quaternion multiplication
%[text] `def == 0`
%[text] for dual quaternion $\\tilde{\\bf q\n}$ and $\\tilde{\\bf p}\n$
%[text] $\\tilde{\\bf q} \\odot \\tilde{\\bf p}= \\left\[\\begin{array}{cc}\n{\\bf q}\_r\\odot & 0\_{4\\times4}     \\\\\n{\\bf q}\_d\\odot & {\\bf q}\_r\\odot    \\\\\n\\end{array}\\right\]$
%[text] where
%[text] ${\\bf q}\\odot{\\bf p}= \\left\[\\begin{array}{cc}\n     q\_{4}\\bar{{\\bf{p}}}+p\_{4}\\bar{{\\bf{q}}} + \\bar{{\\bf{q}}}\\times\\bar{{\\bf{p}}}\\\\\nq\_{4}p\_{4}-\\bar{{\\bf{q}}}^{T}\\bar{{\\bf{p}}}\n\\end{array}\\right\]$
%[text] `def == 1`
%[text] $\\tilde{\\bf q} \\otimes \\tilde{\\bf p}= \\left\[\\begin{array}{cc}\n{\\bf q}\_r\\otimes & 0\_{4\\times4}     \\\\\n{\\bf q}\_d\\otimes & {\\bf q}\_r\\otimes    \\\\\n\\end{array}\\right\]$
%[text] where
%[text] ${\\bf q}\\otimes{\\bf p}= \\left\[\\begin{array}{cc}\n     q\_{4}\\bar{{\\bf{p}}}+p\_{4}\\bar{{\\bf{q}}} - \\bar{{\\bf{q}}}\\times\\bar{{\\bf{p}}}\\\\\nq\_{4}p\_{4}-\\bar{{\\bf{q}}}^{T}\\bar{{\\bf{p}}}\n\\end{array}\\right\]$
%[text] dq: dual quaternion \[real, dual\], nx8
%[text] dp: dual quaternion \[real, dual\], nx8
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20231219  y.yoshimura
%[text] See also qMult.
function dq = dqMult_(scalar, def, dq, dp)

dqReal = dq(:,1:4);
dqDual = dq(:,5:8);

dpReal = dp(:,1:4);
dpDual = dp(:,5:8);

realPart = qMult_(scalar, def, dqReal, dpReal);
dualPart = qMult_(scalar, def, dqReal, dpDual) + qMult_(scalar, def, dqDual, dpReal);

dq = [realPart dualPart];


end


%[appendix]{"version":"1.0"}
%---
