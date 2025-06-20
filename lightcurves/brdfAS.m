%[text] # calculating bidrectional reflectance distribution function (BRDF) of Ashikhmin–Shrley model
%[text] Ashikhmin–Shirleyモデルのdiffuse, specular, NDF計算
%[text] `sat`: satellite configuration read with `readSC`
%[text] `sunB:` sun vector from satellite to Sun expressed with body-fixed frame, Mx3 matrix
%[text] `vB: reference` expressed with body-fixed frame, Mx3 matrix
%[text] cd: diffuse part of BRDF, nFacet x 1 vector
%[text] cs: specular part of BRDF, nFacet x 1 vector
%[text] D: normal distribution function of BRDF, nFacet x 1 vector
%[text] ## note
%[text] Ashikhmin–Shirley model is written as
%[text] $c\_{d} = \\frac{28\\rho}{23 \\pi} \\left(1-F\_{0}\\right) \\left\[1-\\left(1-\\frac{{\\bf{n}}^{T}{\\bf{s}}}{2}\\right)^{5}\\right\]\\left\[1-\\left(1-\\frac{{\\bf{n}}^{T}{\\bf{v}}}{2}\\right)^{5}\\right\]$
%[text] $c\_{s} = \\frac{\\sqrt{\\left(n\_{u}+1\\right)\\left(n\_{v}+1\\right)}}{8\\pi} \\frac{F}{({\\bf{v}}^{T}{\\bf{h}}){\\rm max}({\\bf{n}}^{T}{\\bf{s}},{\\bf{n}}^{T}{\\bf{v}})} D({\\bf{h}})$
%[text] where
%[text] $F = F\_{0} + (1-F\_{0}) (1 - {\\bf{v}}^{T}{\\bf{h}})^{5} \\\\\n$
%[text] $D({\\bf{h}}) = ({\\bf{n}}^{T}{\\bf{h}})^{n\_u\\cos^{2}{\\phi\_{h}} + n\_{v}\\sin^{2}{\\phi\_{h}}}$
%[text] ## references 
%[text] Ashikhmin, Michael, & Shirley, Peter. “An Anisotropic Phong BRDF Model.” Journal of graphics tools, vol. 5, no. 2, 2000, pp. 25-32.
%[text] ## revisions
%[text] 20220317  y.yoshimura, y.yoshimula@gmail.com
%[text] See also srpAS, readSC.
function [cd, cs, D] = brdfAS(sat, sunB, vB) 
sunB = sunB ./ vecnorm(sunB, 2, 2); %一応 normalize, Mx3
vB = vB ./ vecnorm(vB, 2, 2);

h = sunB + vB; % Mx3 matrix, bisector vector of sun and observer vectors
h = h ./ vecnorm(h, 2, 2);

[phiH, thetaH, ~] = cart2sph(h(:,1), h(:,2), h(:,3));
thetaH = pi/2 - thetaH; % +z軸から測った角度にする
phiH = phiH'; % row vector, 1xM
% thetaH = thetaH'; % row vector, 1xM

% sat.normal * sunB'は，faceの数 nFacet x 時間履歴の数 M のmatrixになる
NS = sat.normal * sunB'; % nFacetxM
NV = sat.normal * vB'; % nFacetxM
NH = sat.normal * h';
VH = dot(vB,h,2)'; % 1xM

% nFacet x M matrixへ拡張
nu = repmat(sat.nu, 1, size(sunB,1));
nv = repmat(sat.nv, 1, size(sunB,1));
%[text] ## diffuse, nFacet x M
cd = 28 / 23 .* sat.Cd / pi .* (1 - sat.F0) .* (1 - (1 - NS./2).^5) .* (1 - (1 - NV./2).^5);
%[text] ## specualr
F = sat.F0 + (1 - sat.F0) .* (1 - VH).^5; % NxM
M = sqrt((nu + 1) .* (nv + 1)) ./ 8 ./ pi .* F ./ VH ./ max(NS, NV); % NxM
infIndex = isinf(M);
M(infIndex) = 0; % M計算時にmax関数で0割によるinf発生を除去

%[text] when calculate using $\\left(\\bf{n}^{T}\\bf{h}\\right)^{\\frac{n\_{u}(\\bf{h}^{T}\\bf{u}\_{u})^{2} + n\_{v}(\\bf{h}^{T}\\bf{u}\_{v})^{2}}{1-(\\bf{h}^{T}\\bf{n})^{2}}} $
hu = sat.uu * h'; % nSat x M
hv = sat.uv * h'; 
D = NH.^((nu .* hu.^2 + nv .* hv.^2) ./ (1 - NH.^2));
D(isnan(D))= 1;

%[text] when calculate using$D({\\bf{h}}) = ({\\bf{n}}^{T}{\\bf{h}})^{n\_u\\cos^{2}{\\phi\_{h}} + n\_{v}\\sin^{2}{\\phi\_{h}}}$
% D = NH.^(nu .* cos(phiH).^2 + nv .* sin(phiH).^2);

cs = M .* D; % nFacet x M

% if the facet can be seen or not
cd = cd .* (NS > 0) .* (NV > 0);
cs = cs .* (NS > 0) .* (NV > 0);
D = D .* (NS > 0) .* (NV > 0);

end

%[appendix]{"version":"1.0"}
%---
