%[text] # calculating remaining term $M$ in the Cook–Torrance model
%[text] Cook–TorranceモデルのNDF項以外を計算
%[text] `sat`: satellite configuration read with `readSC`
%[text] `v:` reference vector, nFacet x 3 matrix
%[text] `sunB:` sun vector from satellite to Sun expressed with body-fixed frame, nFacet x 3 matrix
%[text] `M`: remaining term
%[text] $M(v)=\frac{G(v)F(v)}{4}$
%[text] ## note
%[text] NA
%[text] ## references 					
%[text] Analytic Approximation of High-Fidelity Solar Radiation Pressure.
%[text] ## revisions
%[text] 20200915  y.yoshimura, y.yoshimula@gmail.com
%[text] See also srpApproxCT.
function M = ctM2(sat, v, sunB)

nFacet = size(sat.area,1);

% bisector vector, nFacet x 3 matrix
% v has the size of nFacet x 3 matrix
h = v + sunB;
h = h ./ vecnorm(h,2,2);
%[text] ## pre calculation
NS = sat.normal * sunB'; % nFacet x 1
NH = dot(sat.normal, h, 2); % nFacet x 1
NV = dot(sat.normal, v, 2); % nFacet x 1
VH = dot(v, h, 2); % nFacet x 1
%[text] ## specular
nest = (1 + sqrt(sat.F0)) ./ (1 - sqrt(sat.F0)); % nFacet x 1
g = sqrt(nest.^2 + VH.^2 - 1); % nFacetx1

tmp1 =  2 * NH .* NV ./ VH;
tmp2 =  2 * NH .* NS ./ VH;

G = min(1, tmp1);
G = min(G, tmp2); % nFacet x 1

tmp1 = (g - VH).^2 / 2 ./ (g + VH).^2;
tmp2 = (1 + (VH .* (g + VH) - 1).^2 ./ (VH .* (g - VH) + 1).^2);
F = tmp1 .* tmp2;

M = G .* F ./ 4;

end

%[appendix]{"version":"1.0"}
%---
