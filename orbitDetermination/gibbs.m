%[text] # Gibbs method of orbit determination from three position vectors
%[text] `r1, r2, r3: geocentric position vectors at t1, t2, and t3`
%[text]  `order must be t1 < t2 < t3`
%[text] mu: Earth's gravitational constant, unit must be consistent with the position vectros
%[text]         重力定数．位置ベクトルと単位を合わせる．e.g., 位置ベクトルがkmの場合，muは${\\rm km^3/s^2}$にする
%[text] ## note
%[text] `r1, r2, r3, their order must be t1 < t2 < t3`
%[text] `r1, r2, r3`は時系列順に並んでいる必要がある
%[text] ## references 
%[text] Curtis, Howard D. Orbital Mechanics for Engineering Students. Butterworth-Heinemann, 2013, p231-237.
%[text] ## revisions
%[text] 20221110  y.yoshimura, y.yoshimula@gmail.com
%[text] See also orbitConst.
function [v1, v2, v3] = gibbs(r1, r2, r3, mu)
% arguments
%     r1 (1,3) {mustBeNumeric}
%     r2 (1,3) {mustBeNumeric}
%     r3 (1,3) {mustBeNumeric}
%     mu (1,1) {mustBeNumeric}
% end

r1Norm = norm(r1);
r2Norm = norm(r2);
r3Norm = norm(r3);

u1 = r1 ./ r1Norm;

C23 = cross(r2, r3);

% if position vectors are not coplanar
if u1 * C23' > 1e-3
    error('position vectors are not coplanr.')
end

%[text] ## calculate N, D, and S vectors
N = r1Norm .* cross(r2, r3) + r2Norm .* cross(r3, r1) + r3Norm .* cross(r1, r2);% 1x3 vector

D = cross(r1, r2) + cross(r2, r3) + cross(r3, r1);% 1x3 vector

S = r1 .* (r2Norm - r3Norm) + r2 .* (r3Norm - r1Norm) + r3 .* (r1Norm - r2Norm);% 1x3 vector

v1 = sqrt(mu / norm(N) / norm(D)) .* (cross(D, r1)./ r1Norm + S); % 1x3 vector
v2 = sqrt(mu / norm(N) / norm(D)) .* (cross(D, r2)./ r2Norm + S); % 1x3 vector
v3 = sqrt(mu / norm(N) / norm(D)) .* (cross(D, r3)./ r3Norm + S); % 1x3 vector


end

%[appendix]{"version":"1.0"}
%---
