%[text] # numerically calculating solar radiation pressure (SRP) forces and torques with Ashikhmin–Shrley model
%[text] Ashikhmin–ShirleyモデルによるSRP計算
%[text] ## input
%[text] `sat`: satellite configuration read with `readSC`
%[text] `sunB:` sun vector from satellite to Sun expressed with body-fixed frame, 1x3 matrix
%[text] `d`: distance between satellite and Sun, m
%[text] `const`: orbital constants
%[text] `nMC`: number of sampling for Monte-Carlo integration
%[text] ## output
%[text] `sat.force`: SRP force at each facet expressed with body-fixed frame, N, nx3 matrix
%[text] `sat.torque`: SRP torque at each facet expressed with body-fixed frame, Nm, nx3 matrix
%[text] `srpCdOut`: total diffuse part of SRP, N, 1x3
%[text] `srpCsOut`: total specular part of SRP, N, 1x3
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
function [sat, srpCdOut, srpCsOut] = srpAS(sat, sunB, d, const, nMC)
arguments (Input)
    sat
    sunB (:,3) {mustBeNumeric}
    d (:,1) {mustBeNumeric}
    const
    nMC = 10^4
end
arguments (Output)
    sat
    srpCdOut (1,3)
    srpCsOut (1,3)
end
%[text] ## parameters
dAU = km2AU(d ./ 10^3, const); % AU
S0 = const.S0; % Solar constant, W/m^2
c = const.c; % light speed, m/s
coeff = -S0 / c / dAU^2;

sunB = sunB ./ norm(sunB); % normalize, 1x3

NS = sat.normal * sunB'; % nFacetx1
%[text] ## diffuse (analytic)
cd1 = 28/23 .* sat.rho ./ pi .* (1 - sat.F0) .* (1 - (1 - NS / 2).^5);  %nFacetx1
srpCd = [0, 0, sum(cd1,1)* 1573/2688 * pi];
%[text] ## speuclar (numerical)
%[text] $c\_{s} = \\frac{\\sqrt{\\left(n\_{u}+1\\right)\\left(n\_{v}+1\\right)}}{8\\pi} \\frac{F}{({\\bf{v}}^{T}{\\bf{h}}){\\rm max}({\\bf{n}}^{T}{\\bf{s}},{\\bf{n}}^{T}{\\bf{v}})} D({\\bf{h}})$
%[text] $F = F\_{0} + (1-F\_{0}) (1 - {\\bf{v}}^{T}{\\bf{h}})^{5} \\\\\n$
%[text] $D({\\bf{h}}) = ({\\bf{n}}^{T}{\\bf{h}})^{n\_u\\cos^{2}{\\phi\_{h}} + n\_{v}\\sin^{2}{\\phi\_{h}}} = ({\\bf n^T h})^{\\frac{n\_u(h^Tu\_u)^2+n\_v(h^T u\_v)^2}{1-(h^Tn)^2}}$
%[text] #### transform to local frame
%[text] 法線ベクトルが${\\bf n} =\[0, 0, 1\]^T\n$になるlocal frameへ変換して計算する
nFacet = size(sat.normal,1); % # of facet

% sun vector and normal vectors at local frame (normal vector is along z-axis)
sLocal = qRotation(4, repmat(sunB, nFacet,1), sat.qlb); % nFacet x 3
% nLocal = qRotation(4, sat.normal, sat.qlb); % (for debug), nFacet x 3

%[text] ## sampling
%[text] importance sampling for ${\\bf{h}} = \[\\cos\\phi\_h \\sin\\theta\_h, \\sin\\phi\_h \\sin\\theta\_h, \\cos\\theta\_h\]^T$
u1 = rand(1,nMC); % 1 x nMC vectors
u2 = rand(1,nMC);

% phiH = atan(sqrt((sat.nu + 1) / (sat.nv + 1)) * tan(2*pi .* u1));
% phiH = phiH .* (u1 <= 0.25) + (phiH + pi) .* (u1 > 0.25) .* (u1 < 0.75)...
%     + (phiH + 2 * pi) .* (u1 >= 0.75);

phiFun = @(u, nu, nv) atan(sqrt((nu + 1) ./ (nv + 1)) .* tan(pi .* u / 2)); % nFacet x nMC
ind1 = (u1 < 0.25);
ind2 = (u1 >= 0.25) .* (u1 < 0.5);
ind3 = (u1 >= 0.5) .* (u1 < 0.75);
ind4 = (u1 >= 0.75);

u1 = ind1 .* 4 .* u1 +  ind2 .* (4 * (u1 -0.25)) ...
    +  ind3 .* (4 * (u1 - 0.5)) + ind4 .* (4 * (u1 - 0.75));

phiH = phiFun(u1, sat.nu, sat.nv);
phiH = ind1 .* phiH  + ind2 .* (phiH + pi/2) ...
    + ind3 .* (phiH + pi) + ind4 .* (phiH + 3/2*pi); %nFacet x nMC

thetaH = acos(u2.^(1 ./ (sat.nu .* cos(phiH).^2 + sat.nv .* sin(phiH).^2 + 1)));  %nFacet x nMC

% h = [sin(thetaH).*cos(phiH) sin(thetaH).*sin(phiH) cos(thetaH)]; %Nx3

hx = sin(thetaH).*cos(phiH); %nFacet x nMC
hy = sin(thetaH).*sin(phiH);
hz = cos(thetaH);

%[text] ### integration        
% vx = 2 * (sLocal * h') .* h(:,1)' - sLocal(:,1); % nFacet x nMC
% vy = 2 * (sLocal * h') .* h(:,2)' - sLocal(:,2); % nFacet x nMC
% vz = 2 * (sLocal * h') .* h(:,3)' - sLocal(:,3); % nFacet x nMC

SH = sLocal(:,1) .* hx + sLocal(:,2) .* hy + sLocal(:,3) .* hz; % nFacet x nMC

vx = 2 * SH .* hx - sLocal(:,1); % nFacet x nMC
vy = 2 * SH .* hy - sLocal(:,2);
vz = 2 * SH .* hz - sLocal(:,3);

VH = vx .* hx + vy .* hy + vz .* hz; % nFacet x nMC
NV = vz; % nFacet x nMC
NH = hz; % nFacet x nMC

F = sat.F0 + (1 - sat.F0) .* (1 - VH).^5; % nFacet x nMC

M = 1 ./ VH ./ max(NS, NV); % nFacet x nMC
M(isinf(M)) = 0; % Inf項を消す

% weight
W = abs(SH) .* NV ./ NH .* M;

temp = W .* F; % nxnMC

% for SRP
csASx = (NV > 0) .* temp .* vx; % nFacetxN matrix
csASy = (NV > 0) .* temp .* vy; % nFacetxN matrix
csASz = (NV > 0) .* temp .* vz; % nFacetxN matrix

% specular component, nx3
tmp = [mean(csASx,2), mean(csASy,2), mean(csASz,2)];
srpCs = qRotation(4, tmp, qInv(4, sat.qlb)); % transform back to the original body-fixed frame, nFacet x 3

%[text] ## shadowing
%[text] 1: sunlit, 0: shade
sunlitFlag = (NS > 0); % nFacetx1 matrix, 1: sunlit, 0: shade

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
