%[text] # calculating remaining term $M$ in the Cook–Torrance model
%[text] Cook–TorranceモデルのNDF項以外を計算
%[text] `sat`: satellite configuration read with `readSC`
%[text] `v:` reference vector, 1x3 or 3x1 vector
%[text] `sunB:` sun vector from satellite to Sun expressed with body-fixed frame, 1x3 or 3x1 vector
%[text] `M`: remaining term
%[text] $M({\\bf v})=\\frac{G({\\bf v})F({\\bf v})}{4} $
%[text] ## note
%[text] `sat.normal` may be nFacet x 3 and `sat.F0` may be scalar or nFacet x 1.
%[text] ## references
%[text] Analytic Approximation of High-Fidelity Solar Radiation Pressure.
%[text] ## revisions
%[text] 20200915  y.yoshimura, y.yoshimula@gmail.com
%[text] 20260706  consolidated duplicate implementations (srp/ctM.m, hifiSRP/ctM.m)
%[text] See also srpApproxCT.
function M = ctM(sat, v, sunB)

n = sat.normal; % normal vectors, nx3 matrix

v = v(:); % 3x1 vector
sunB = sunB(:); % 3x1 vector

% bisector vector, 3x1 vector
h = v + sunB;
h = h ./ norm(h);
%[text] ## specular
nest = (1 + sqrt(sat.F0)) ./ (1 - sqrt(sat.F0));
g = sqrt(nest.^2 + (v' * h).^2 - 1); % nx1

temp1 =  2 * (n * h) .* (n * v) ./ (v' * h);
temp2 =  2 * (n * h) .* (n * sunB) ./ (v' * h); % nx1

G = min(1, temp1);
G = min(G, temp2); % nx1

temp1 = (g - v' * h).^2 / 2 ./ (g + v' * h).^2;
temp2 = (1 + ((v' * h) .* (g + v' * h) - 1).^2 ./ ((v' * h) .* (g - v' * h) + 1).^2);
F = temp1 .* temp2;

M = sum(G .* F, 1) ./ 4;

end

%[appendix]{"version":"1.0"}
%---
