%[text] # calculating optimal continuous quaternion from direction cosine matrix (DCM)
%[text] `scalar:` specify the definition of the quaternion 
%[text] `scalar == 0`
%[text] ${\\bf q} = \[q\_0,q\_1,q\_2,q\_3\]^T=\[\\cos(\\theta/2), {\\bf e}^T\\sin(\\theta/2)\]^T$
%[text] `scalar == 4`
%[text] ${\\bf q} = \[q\_1,q\_2,q\_3, q\_4\]^T=\[{\\bf e}^T\\sin(\\theta/2), \\cos(\\theta/2)\]^T$
%[text] `qK`: previous quaternion, $q\_{k-1}\n$1x4 matrix
%[text] dcm: direction cosine matrix, 3x3 matrix
%[text] `q`: quaternion, 1x4 matrix
%[text] ## note
%[text] NA
%[text] ## references 
%[text] Wu, Jin. “Optimal Continuous Unit Quaternions from Rotation Matrices.” Journal of Guidance, Control, and Dynamics, vol. 42, no. 4, 2019, pp. 919-22.
%[text] ## revisions
%[text] 20221008  y.yoshimura, y.yoshimula@gmail.com
%[text] See also dcm2q.
function q = dcm2qContinuous(scalar, qK, dcm)
% arguments
%     scalar (1,1) {mustBeMember(scalar, [0, 4])}
%     qK (1,4) {mustBeNumeric}
%     dcm (3,3) {mustBeNumeric}
% end

% temporarily the definition that q4 = cos(theta/2) is used
if scalar == 0
    qK = [qK(:,2), qK(:,3), qK(:,4), qK(:,1)];
end

% Eq. (7)
B = 1/3 .* dcm;
z = [ B(2,3) - B(3,2)
    B(3,1) - B(1,3)
    B(1,2) - B(2,1)];
K = [B + B' - trace(B) .* eye(3) z
    z' trace(B)];
[V, ~] = eig(K);

q = V * diag([0 0 0 1]) * V^(-1) * qK';

q = q' ./ norm(q);

if scalar == 0
    q = [q(4) q(1) q(2) q(3)];
end

end

%[appendix]{"version":"1.0"}
%---
