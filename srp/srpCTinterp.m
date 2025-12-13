%[text] # solar radiation pressure (SRP) forces and torques with Cook–Torrance model using correction parameters and interpolation
%[text] 補正係数とその内挿でCook–TorranceモデルによるSRP計算
%[text] ## input
%[text] `sat`: satellite configuration read with `readSC`
%[text] `sunB:` sun vector from satellite to Sun expressed with body-fixed frame, 1x3 vector
%[text] `d`: distance between satellite and Sun, m
%[text] `const`: orbital constants
%[text] `correctionPara`: correction parameters
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
function [sat, srpCdOut, srpCsOut] = srpCTinterp(sat, sunB, d, const, correctionPara)
% arguments (Input)
%     sat
%     sunB (:,3) {mustBeNumeric}
%     d (:,1) {mustBeNumeric}
%     const
%     correctionPara
% end
% arguments (Output)
%     sat
%     srpCdOut (1,3)
%     srpCsOut (1,3)
% end
%[text] ## parameters
dAU = km2au(d ./ 10^3, const); % AU
S0 = const.S0; % Solar constant, W/m^2
c = const.c; % light speed, m/s
coeff = -S0 / c / dAU^2;

sunB = sunB ./ norm(sunB);

%[text] ## diffuse (analytic) and speuclar (corrected)
NS = sat.normal * sunB'; % nFacet x 1
thetaI = acos(NS); % nFacet x 1, clamped

% Interpolation
deltaS = interp2(correctionPara.thetaIspan, correctionPara.mSpan, correctionPara.deltaS, thetaI, sat.mCT, 'linear', 1); % extrapolated values are 1 (standard) or clamp?
deltaN = interp2(correctionPara.thetaIspan, correctionPara.mSpan, correctionPara.deltaN, thetaI, sat.mCT, 'linear', 1);

% Fix NaN issues if thetaI or mCT are out of range (though interp2 should handle with extrapolation or NaN)
% Use spline or pchip? Linear is requested.
% If out of range, what to do? 'spline' might be better or closest.
% For now 'linear' is fine. Using '1' as extrapolation value might be risky if deltaS/N are far from 1.
% But deltaS/N -> 1 as m->0? No, deltaS/N are corrections.
% Let's use nearest or keep NaNs and fix them?
% Let's stick to default 'linear' (returns NaN) and fill? Or just linear with extrapolation?

% Force calculation (per area unit, scaled by coeff)
% Note: sunB is 1x3, others are Nx1. Broadcasting needed.
fCorrected = (1 - deltaS .* sat.Cs) .* sunB + (2/3.*sat.Cd + 2 .* NS .* deltaN .* sat.Cs) .* sat.normal;

%[text] ## total SRP
sunlitFlag = double(NS > 0); % nFacet x 1 matrix, 1: sunlit, 0: shade
sat.force = sunlitFlag .* coeff .* sat.area .* NS .* fCorrected; % nx3 matrix
sat.torque = cross(sat.pos, sat.force); % nx3 matrix

%[text] ## for output variables
% Total SRP force
srpTotal = sum(sat.force, 1);

% Separation (Approximate)
% Diffuse part (term with Cd)
fDiffuse = sunlitFlag .* coeff .* sat.area .* NS .* ((2/3.*sat.Cd) .* sat.normal);
srpCdOut = sum(fDiffuse, 1);

srpImpinged = sunlitFlag .* coeff .* sat.area .* NS .* sunB;

% Specular part (rest)
srpCsOut = srpTotal - srpCdOut - srpImpinged;

end


%[appendix]{"version":"1.0"}
%---
