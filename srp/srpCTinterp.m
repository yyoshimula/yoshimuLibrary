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
function [sat, srpCdOut, srpCsOut] = srpCTinterp(sat, sunB, d, const, correctionPara, mCase)
% arguments (Input)
%     sat
%     sunB (:,3) {mustBeNumeric}
%     d (:,1) {mustBeNumeric}
%     const
%     NDF = 'Beckmann' % 現状Beckamnn distributionのみ対応
%     nMC = 10^3;
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
thetaI = acos(sat.normal * sunB'); % nFacet x 1
fSRP = (1 - deltaS(mCase,thetaI) * sat.Cs) .* sunB + (2/3*sat.Cd + 2 * NS * deltaN(i,j) * sat.Cs) .* sat.normal;
fCorrected(j,:,i) = coeff * fCorrected(j,:,i);
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
