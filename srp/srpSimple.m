%[text] # calculating solar radiation pressure (SRP) forces and torques with Lambertian and perfect specular model
%[text] 完全ランバート反射，完全鏡面反射モデルによるSRP計算
%[text] ## input
%[text] `sat`: satellite configuration read with `readSC`
%[text] `sunB:` sun vector from satellite to Sun expressed with body-fixed frame, 1x3 unit vector
%[text] `d`: distance between satellite and Sun, m
%[text] `const`: orbital constants
%[text] ## output
%[text] `sat.force`: SRP force at each facet expressed with body-fixed frame, N, nx3 matrix
%[text] `sat.torque`: SRP torque at each facet expressed with body-fixed frame, Nm, nx3 matrix
%[text] `srpCdOut`: total diffuse SRP force expressed with body-fixed frame, N, 1x3 matrix
%[text] `srpCsOut`: total specular SRP force expressed with body-fixed frame, Nm, 1x3 matrix
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] See also readSC, orbitConst.
function [sat, srpCdOut, srpCsOut] = srpSimple(sat, sunB, d, const) %#codegen
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
S0 = const.S0; % Solar constant, W/m^2
c = const.c; % light speed, m/s

Bf = 2/3;
kappa = sat.kappa; % thermal emissivityを考慮するときは変更

dAU = km2AU(d ./ 10^3, const); % AU
sunB = sunB(:); % column vector
sunB = sunB ./ norm(sunB);

coeff = -S0 / c / dAU^2;
%[text] ## shadowing flag
sunlitFlag = (sat.normal * sunB > 0); % nx1 matrix, 1: sunlit, 0: shade
ns = sat.normal * sunB; % nx1 matrix
%[text] ## calculate
%[text] satellite has $n$ facets
%[text] ### diffuse
srpCd = 2 / 3 * sat.Cd .* sat.normal; % nx3
%[text] ### specular
rRef = 2 * ns .* sat.normal - sunB'; % nx3
srpCs = sat.Cs .* rRef;
%[text] ### total
tmp = ns .* (sat.Ca + sat.Cd) .* sunB' ...
    + ns .* (Bf * sat.Cd + kappa .* sat.Ca ...
    + 2.0 .* sat.Cs .* ns) .* sat.normal; % nx3 matrix
sat.force = sunlitFlag .* coeff .* sat.area .* tmp; % nx3 matrix
sat.torque = cross(sat.pos, sat.force); % nx3 matrix
%[text] ## for output variables
%[text] diffuse part of SRP and specular part of SRP
tmp = coeff .* sat.area .* ns .* sunlitFlag .* srpCd; % nFacet x 3
srpCdOut = sum(tmp,1); % 1x3

tmp = coeff .* sat.area .* ns .* sunlitFlag .* srpCs;
srpCsOut = sum(tmp,1); % 1x3

end

%[appendix]{"version":"1.0"}
%---
