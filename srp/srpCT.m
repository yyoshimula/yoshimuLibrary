%[text] # numerically calculating solar radiation pressure (SRP) forces and torques with Cook–Torrance model using importance sampling
%[text] Cook–TorranceモデルによるSRP計算
%[text] ## input
%[text] `sat`: satellite configuration read with `readSC`
%[text] `sunB:` sun vector from satellite to Sun expressed with body-fixed frame, 1x3 vector
%[text] `d`: distance between satellite and Sun, m
%[text] `const`: orbital constants
%[text] NDF: (optional) NDF distribution function, default: Beckamnn distribution
%[text] `nMC`: (optional) number of random numbers for integration
%[text] ## output
%[text] `sat.force`: SRP force at each facet expressed with body-fixed frame, N, nx3 matrix
%[text] `sat.torque`: SRP torque at each facet expressed with body-fixed frame, Nm, nx3 matrix
%[text] `srpCdOut`: total diffuse part of SRP, N, 1x3 vector
%[text] `srpCsOut`: total specular part of SRP, N, 1x3 vector
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] See also srpLps, readSC, orbitConst.
function [sat, srpCdOut, srpCsOut] = srpCT(sat, sunB, d, const, NDF, nMC)
arguments (Input)
    sat
    sunB (:,3) {mustBeNumeric}
    d (:,1) {mustBeNumeric}
    const
    NDF = 'Beckmann' % 現状Beckamnn distributionのみ対応
    nMC = 10^3;
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
%[text] ## diffuse (analytic)
srpCd = 2 / 3 * sat.rho .* sat.normal; % nx3
%[text] ## speuclar (numerical)
%[text] #### transform to local frame
%[text] 法線ベクトルが${\\bf n} =\[0, 0, 1\]^T\n$になるlocal frameへ変換して計算する
nFacet = size(sat.normal,1); % # of facet

% sun vector and normal vectors at local frame (normal vector is along z-axis)
sLocal = qRotation(4, repmat(sunB,nFacet,1), sat.qlb); % nFacet x 3
% nLocal = qRotation(4, sat.normal, qlb); % (for debug), nFacet x 3

NS = sat.normal * sunB'; % nx1

%[text] ## Monte-Carlo integration
if strcmp(NDF, 'Beckmann')
%[text] ## sampling for Beckmann distribution
%[text] importance sampling for ${\\bf{h}} = \[\\cos\\phi\_h \\sin\\theta\_h, \\sin\\phi\_h \\sin\\theta\_h, \\cos\\theta\_h\]^T$
    u1 = rand(nMC,1);
    thetaH = atan(sqrt(-sat.mCT.^2 .* log(u1'))); % nFacet x N
    phiH = 2 * pi .* rand(1, nMC); % 1xN

    hx = sin(thetaH).*cos(phiH);
    hy = sin(thetaH).*sin(phiH);
    hz = cos(thetaH); % half vectors, nFacet x N
else % for Gauss distribution
    phiR = 2 * pi .* rand(1, nMC);
    thetaR = pi / 2 .* rand(1, nMC);

    vx = sin(thetaR).*cos(phiR); % 1xnMC
    vy = sin(thetaR).*sin(phiR);
    vz = cos(thetaR);

    hx = sLocal(:,1) + vx;
    hy = sLocal(:,2) + vy;
    hz = sLocal(:,3) + vz; % nFacet x nMC

    % Normalize the half vectors
    hx = hx ./ sqrt(hx.^2 + hy.^2 + hz.^2);
    hy = hy ./ sqrt(hx.^2 + hy.^2 + hz.^2);
    hz = hz ./ sqrt(hx.^2 + hy.^2 + hz.^2);

    thetaH = acos(hz);

    D = exp(-(thetaH ./ sat.mCT).^2); % Gaussian distribution, nxN
end
%[text] ### integration
SH = sLocal(:,1) .* hx + sLocal(:,2) .* hy + sLocal(:,3) .* hz;% nFacet x N
vx = 2 * SH .* hx - sLocal(:,1); % nFacet x N
vy = 2 * SH .* hy - sLocal(:,2); % nFacet x N
vz = 2 * SH .* hz - sLocal(:,3); % nFacet x N

VH = vx .* hx + vy .* hy + vz .* hz; % nFacet x N
NV = vz; % nFacet x N
NH = hz; % nFacet x N

nest = (1 + sqrt(sat.F0)) ./ (1 - sqrt(sat.F0)); % nFacet x 1
g = sqrt(nest.^2 + VH.^2 - 1); % nFacet x N

temp1 = 2 * NH .* NV ./ VH; % nFacet x N
temp2 = 2 * NH .* NS ./ VH; % nFacet x N
G = min(1, temp1);
G = min(G, temp2); % nFacet x N

temp1 = (g - VH).^2 / 2 ./ (g + VH).^2;
temp2 = (1 + (VH .* (g + VH) - 1).^2 ./ (VH .* (g - VH) + 1).^2);
F = temp1 .* temp2; % nFacet x N

% weight
W = abs(SH) .* G ./ NS ./ NH; % nFacet x N
tmp = W .* F; % nFacet x N
tmp(isnan(tmp)) = 0; % NaNを0にする

% SRP
csCTx = (NV > 0) .* tmp .* vx; % nFacet x N  matrix
csCTy = (NV > 0) .* tmp .* vy; % nFacet x N  matrix
csCTz = (NV > 0) .* tmp .* vz; % nFacet x N  matrix

% specular component, nx3
tmp = [mean(csCTx, 2), mean(csCTy, 2), mean(csCTz, 2)];

srpCs = qRotation(4, tmp, qInv(4, sat.qlb)); % transform back to the original body-fixed frame, nFacet x 3

%[text] ## total SRP
sunlitFlag = (NS > 0); % nFacet x 1 matrix, 1: sunlit, 0: shade
tmp = (sunB + srpCd + srpCs);
sat.force = sunlitFlag .* coeff .* sat.area .* NS .* tmp; % nx3 matrix
sat.torque = cross(sat.pos, sat.force); % nx3 matrix
%[text] ## for output variables
%[text] diffuse part of SRP and specular part of SRP
% Calculate the contribution of each facet to the solar radiation pressure coefficient
tmp = coeff .* sat.area .* NS .* sunlitFlag .* srpCd; % nFacet x 3
% Sum contributions across facets to get total for srpCd
srpCdOut = sum(tmp,1); % 1x3

% Calculate the contribution of each facet to the solar radiation pressure coefficient for srpCs
tmp = coeff .* sat.area .* NS .* sunlitFlag .* srpCs;
% Sum contributions across facets to get total for srpCs
srpCsOut = sum(tmp,1);
end

%[appendix]{"version":"1.0"}
%---
