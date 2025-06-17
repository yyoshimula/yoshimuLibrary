%[text] # numerically calculating solar radiation pressure (SRP) forces and torques with Cook–Torrance model using uniform distribution (not recommended)
%[text] Cook–TorranceモデルによるSRP計算
%[text] `sat`: satellite configuration read with `readSC`
%[text] `sunB:` sun vector from satellite to Sun expressed with body-fixed frame, 1x3 vector
%[text] `d`: distance between satellite and Sun, m
%[text] `const`: orbital constants
%[text] `nMC`: number of random numbers for integration
%[text] `varagin`: NDF distribution function, default: Beckamnn distribution, or Gauss distribution if specified
%[text] `sat.force`: SRP force at each facet expressed with body-fixed frame, N, nx3 matrix
%[text] `sat.torque`: SRP torque at each facet expressed with body-fixed frame, Nm, nx3 matrix
%[text] `srpCdOut`: diffuse part of SRP
%[text] `srpCsOut`: specular part of SRP
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20220201  y.yoshimura, y.yoshimula@gmail.com
%[text] See also srpLps, readSC, orbitConst.
function [sat, srpCdOut, srpCsOut] = srpCTuni(sat, sunB, d, const, NDF, nMC)
arguments (Input)
    sat
    sunB (:,3) {mustBeNumeric}
    d (:,1) {mustBeNumeric}
    const
    NDF = 'Beckmann' % Gaussという文字列を入れた場合Gaussian NDFを使う
    nMC = 10^6;
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

sunB = sunB ./ norm(sunB);
%[text] ## sampling
%[text] uniform sampling for spherical integration （**not hemispherical**)
%[text] ${\\bf{v}} = \[\\cos\\phi\_r \\sin\\theta\_r, \\sin\\phi\_r \\sin\\theta\_r, \\cos\\theta\_r\]^T$
% Monte-Carlo integration condition

% sampling for 全球
thetaR = pi .* rand(nMC,1); % nMC x 1
phiR = 2 * pi .* rand(nMC,1);

v = [sin(thetaR).*cos(phiR) sin(thetaR).*sin(phiR) cos(thetaR)]; % nMCx3
h = repmat(sunB, nMC, 1) + v;
h = h ./ vecnorm(h,2,2); % bisector vector, nMCx3
thetaH = acos(sat.normal * h')'; % thetaH, nMCx1

VH = dot(v, h, 2)'; % 1xnMC
NV = sat.normal * v'; % nFacet x nMC, 列方向にreference vector, vの値
NS = sat.normal * sunB'; % nFacet x 1
NH = sat.normal * h'; % nFacet x nMC
%[text] ## diffuse (analytic)
srpCd = 2 / 3 * sat.rho .* sat.normal; % nFacet x 3
%[text] ## speuclar (numerical)
if strcmp(NDF, 'Beckmann')
    D = exp(-(tan(thetaH')./sat.mCT).^2); % Beckmann distribution, nFacet x nMC
    D = D ./ pi ./ sat.mCT.^2 ./ cos(thetaH').^4;
else
    D = exp(-(thetaH' ./ sat.mCT).^2); % Gaussian distribution, nFacet x nMC
end
nest = (1 + sqrt(sat.F0)) ./ (1 - sqrt(sat.F0)); % nx1
g = sqrt(nest.^2 + VH.^2 - 1); % nxN

temp1 = 2 * NH .* NV ./ VH; % nxN
temp2 = 2 * NH .* NS ./ VH; % nxN
G = min(1, temp1);
G = min(G, temp2); % nxN

temp1 = (g - VH).^2 / 2 ./ (g + VH).^2;
temp2 = (1 + (VH .* (g + VH) - 1).^2 ./ (VH .* (g - VH) + 1).^2);
F = temp1 .* temp2; % nxN

temp = D .* G .* F ./ NS ./ NV / 4; % nxN

temp = temp .* NV .* sin(thetaR)'; % nxN
csCTx = (NV > 0) .* temp .* v(:,1)'; % nxN matrix
csCTy = (NV > 0) .* temp .* v(:,2)'; % nxN matrix
csCTz = (NV > 0) .* temp .* v(:,3)'; % nxN matrix

% probability
pCT = 1 / pi * 1 / (2 * pi); % scalar

% specular component, nx3
srpCs = [sum(csCTx, 2), sum(csCTy, 2), sum(csCTz, 2)];
srpCs = srpCs ./ pCT ./ nMC;
%[text] ## shadowing
%[text] 1: sunlit, 0: shade
sunlitFlag = (NS > 0); % nx1 matrix, 1: sunlit, 0: shade

% satellite has n facets
temp = (sunB + srpCd + srpCs);
sat.force = sunlitFlag .* coeff .* sat.area .* NS .* temp; % nx3 matrix
sat.torque = cross(sat.pos, sat.force); % nx3 matrix
%[text] ## for output variables
%[text] diffuse part of SRP and specular part of SRP
tmp = coeff .* sat.area .* NS .* sunlitFlag .* srpCd; % nFacet x 3
srpCdOut = sum(tmp,1); % 1x3

tmp = coeff .* sat.area .* NS .* sunlitFlag .* srpCs;
srpCsOut = sum(tmp,1);
end

%[appendix]{"version":"1.0"}
%---
