%[text] # converting generalized Rodrigues parameters to quaternions 
%[text] `scalar == 0`
%[text] ${\\bf q} = \[q\_0,q\_1,q\_2,q\_3\]^T=\[\\cos(\\theta/2), {\\bf e}^T\\sin(\\theta/2)\]^T$
%[text] `scalar == 4`
%[text] ${\\bf q} = \[q\_1,q\_2,q\_3, q\_4\]^T=\[{\\bf e}^T\\sin(\\theta/2), \\cos(\\theta/2)\]^T$
%[text] `f:` scaling parameter for Rodrigues parameters, scalar
%[text] `a:` parameter for Rodrigues parameters, scalar 
%[text] `p`, generalized Rodrigues parameters, nx3 matrix
%[text] `q`, quaternions, nx4 matrix
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20190121  y.yoshimura, y.yoshimula@gmail.com
%[text] See also GRP2q.
function q = grp2q(scalar, f, a, p)
% arguments
%     scalar (1,1) {mustBeMember(scalar, [0, 4])}
%     f (1,1) {mustBeNumeric}
%     a (1,1) {mustBeNumeric}
%     p (:,3) {mustBeNumeric}
% end

q4Num = -a .* vecnorm(p, 2, 2).^2 + f .* sqrt(f^2 + (1-a^2) .* vecnorm(p, 2, 2).^2);
q4Den = f^2 + vecnorm(p, 2, 2).^2;

q4 = q4Num ./ q4Den;

qv = (a + q4) .* p ./ f;

if scalar == 0
    q = [q4, qv];
else
    q = [qv, q4];
end

end

%[appendix]{"version":"1.0"}
%---
