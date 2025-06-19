%[text] # quaternion averaging with scalar weights
%[text] the definition of output quaternion is consistent with the input quaternion
%[text] 入力のquaternionと同じ定義のquaternionが出力される
%[text] ${\\bf q} = \[q\_0,q\_1,q\_2,q\_3\]^T=\[\\cos(\\theta/2), {\\bf e}^T\\sin(\\theta/2)\]^T$
%[text] ${\\bf q} = \[q\_1,q\_2,q\_3, q\_4\]^T=\[{\\bf e}^T\\sin(\\theta/2), \\cos(\\theta/2)\]^T$
%[text] `q`: quaternions, nx4 matrix
%[text] w: scalar weights
%[text] ## note
%[text] NA
%[text] ## references 
%[text] Markley, F. L., Cheng, Y., Crassidis, J. L., & Oshman, Y. (2007). Averaging Quaternions. Journal of Guidance, Control, and Dynamics, 30(4), 1193-1197. https://doi.org/10.2514/1.28949
%[text] ## revisions
%[text] 20230614  y.yoshimura, y.yoshimula@gmail.com
%[text] See also q2zyx.
function qAveraged = qAve(q, w)
% arguments    
%     q (:,4) {mustBeNumeric}
%     w (:,1) {mustBeNumeric} = ones(size(q,1), 1) ./ size(q,1)
% end

M = zeros(4,4);

for i = 1:size(q,1)
    M = M + w(i) .* (q(i,:)' * q(i,:));
end

[u, ~] = eig(M);

qAveraged = u(:,4)';

end

%[appendix]{"version":"1.0"}
%---
