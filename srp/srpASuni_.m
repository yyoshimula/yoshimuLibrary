%[text] # numerically calculating solar radiation pressure (SRP) forces and torques with Ashikhmin–Shrley model (uniform sampling, not recommended)
%[text] Ashikhmin–ShirleyモデルによるSRP計算
%[text] `sat`: satellite configuration read with `readSC`
%[text] `sunB:` sun vector from satellite to Sun expressed with body-fixed frame, 1x3 matrix
%[text] `d`: distance between satellite and Sun, m
%[text] `const`: orbital constants
%[text] `nMC`: number of sampling for Monte-Carlo integration
%[text] `sat.force`: SRP force at each facet expressed with body-fixed frame, N, nx3 matrix
%[text] `sat.torque`: SRP torque at each facet expressed with body-fixed frame, Nm, nx3 matrix
%[text] `srpCdOut`: diffuse part of SRP
%[text] `srpCsOut`: specular part of SRP
%[text] ## note
%[text] SRP is written as
%[text] $ -\\frac{S\_{0}}{cr^{2}\_{\\rm AU}}A({\\bf{n}}^{T}{\\bf{s}})\\left({{\\bf{s}}}+\\int \\int f\_{r}({\\bf{n}}^{T}{\\bf{v}})\\sin{\\theta\_{r}}{{\\bf{v}}}\\mathrm{d}\\theta\_{r}\\mathrm{d}\\phi\_{r}\\right)$
%[text] Ashikhmin–Shirley model is written as
%[text] $c\_{d} = \\frac{28\\rho}{23 \\pi} \\left(1-F\_{0}\\right) \\left\[1-\\left(1-\\frac{{\\bf{n}}^{T}{\\bf{s}}}{2}\\right)^{5}\\right\]\\left\[1-\\left(1-\\frac{{\\bf{n}}^{T}{\\bf{v}}}{2}\\right)^{5}\\right\]$
%[text] $c\_{s} = \\frac{\\sqrt{\\left(n\_{u}+1\\right)\\left(n\_{v}+1\\right)}}{8\\pi} \\frac{F}{({\\bf{v}}^{T}{\\bf{h}}){\\rm max}({\\bf{n}}^{T}{\\bf{s}},{\\bf{n}}^{T}{\\bf{v}})} D({\\bf{h}})$
%[text] where
%[text] $F = F\_{0} + (1-F\_{0}) (1 - {\\bf{v}}^{T}{\\bf{h}})^{5} \\\\\n$
%[text] $D({\\bf{h}}) = ({\\bf{n}}^{T}{\\bf{h}})^{n\_u\\cos^{2}{\\phi\_{h}} + n\_{v}\\sin^{2}{\\phi\_{h}}}$
%[text] ## references 
%[text] Ashikhmin, Michael, & Shirley, Peter. “An Anisotropic Phong BRDF Model.” Journal of graphics tools, vol. 5, no. 2, 2000, pp. 25-32.
%[text] See also srpLps, readSC, orbitConst.
function [sat, srpCdOut, srpCsOut] = srpASuni_(sat, sunB, d, const, nMC)
%[text] ## parameters
dAU = km2AU(d ./ 10^3, const); % AU
S0 = const.S0; % Solar constant, W/m^2
c = const.c; % light speed, m/s
coeff = -S0 / c / dAU^2;

sunB = sunB ./ norm(sunB); % normalize, 1x3

%% condition
% nFacet = size(sat.normal,1);

%[text] ## sampling
%[text] uniform sampling for spherical integration （**not hemispherical**)
%[text] ${\\bf{v}} = \[\\cos\\phi\_r \\sin\\theta\_r, \\sin\\phi\_r \\sin\\theta\_r, \\cos\\theta\_r\]^T$
thetaR = pi .* rand(nMC,1);
phiR = 2 * pi .* rand(nMC,1);

v = [sin(thetaR).*cos(phiR) sin(thetaR).*sin(phiR) cos(thetaR)]; % nMCx3
h = sunB + v; % bisector
h = h ./ vecnorm(h,2,2); % nMCx3

NH = sat.normal * h'; % nFacet x nMC
NS = sat.normal * sunB'; % nFacetx1
VH = dot(v, h, 2)'; % 1xnMC
NV = sat.normal * v'; % nFacet x nMC, 列方向に各reference vector, vの値
%[text] ## diffuse (analytic)
cd1 = 28/23 .* sat.rho ./ pi .* (1 - sat.F0) .* (1 - (1 - NS / 2).^5);  %nFacetx1
srpCd = [0, 0, sum(cd1,1)* 1573/2688 * pi];
%[text] ## speuclar (numerical)
%[text] $c\_{s} = \\frac{\\sqrt{\\left(n\_{u}+1\\right)\\left(n\_{v}+1\\right)}}{8\\pi} \\frac{F}{({\\bf{v}}^{T}{\\bf{h}}){\\rm max}({\\bf{n}}^{T}{\\bf{s}},{\\bf{n}}^{T}{\\bf{v}})} D({\\bf{h}})$
%[text] $F = F\_{0} + (1-F\_{0}) (1 - {\\bf{v}}^{T}{\\bf{h}})^{5} \\\\\n$
%[text] $D({\\bf{h}}) = ({\\bf{n}}^{T}{\\bf{h}})^{n\_u\\cos^{2}{\\phi\_{h}} + n\_{v}\\sin^{2}{\\phi\_{h}}} = ({\\bf n^T h})^{\\frac{n\_u(h^Tu\_u)^2+n\_v(h^T u\_v)^2}{1-(h^Tn)^2}}$
F = sat.F0 + (1 - sat.F0) .* (1 - VH).^5; % nFacet x nMC
k1 = sqrt((sat.nu+1) .* (sat.nv+1)) / 8 / pi; % nFacet x 1
M = 1 ./ VH ./ max(NS, NV); % nFacet x nMC
M(isinf(M)) = 0; % Inf項を消す

k2 = (sat.nu .* (sat.uu * h').^2 + sat.nv .* (sat.uv * h').^2) ./ (1 - NH.^2);
D = NH .^ k2; % nFacet x nMC

tmp = k1 .* F .* M .* D .* NV .* sin(thetaR)'; % nFacet x nMC

csASx = (NV > 0) .* tmp .* v(:,1)'; % nFacet x nMC matrix
csASy = (NV > 0) .* tmp .* v(:,2)';
csASz = (NV > 0) .* tmp .* v(:,3)';

pdfAS = 1 / pi * 1 / (2 * pi); % probability

srpCs = [sum(csASx, 2), sum(csASy, 2), sum(csASz, 2)]; %nFacet x 3
srpCs = srpCs ./ pdfAS ./ nMC; % nFacet x 3
%[text] ## shadowing
%[text] 1: sunlit, 0: shade
sunlitFlag = (sat.normal * sunB' > 0); % nFacetx1 matrix, 1: sunlit, 0: shade

tmp = sunlitFlag .* (sunB + srpCd + srpCs); % nFacet x 3
sat.force = coeff .* sat.area .* NS .* tmp; % nFacet x 3 matrix
sat.torque = cross(sat.pos, sat.force); % nFacet x 3 matrix
%[text] ## for output variables
%[text] diffuse part of SRP and specular part of SRP
tmp = coeff .* sat.area .* NS .* sunlitFlag .* srpCd; % nFacet x 3
srpCdOut = sum(tmp,1); % 1x3

tmp = coeff .* sat.area .* NS .* sunlitFlag .* srpCs;
srpCsOut = sum(tmp,1);
end

%[appendix]{"version":"1.0"}
%---
