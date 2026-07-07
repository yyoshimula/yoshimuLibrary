%[text] # approximating specular term of SRP with Cook–Torrance model
%[text] Cook–Torranceモデルを用いたSRP近似解析解
%[text] `sat`: satellite configuration read with `readSC`
%[text] `thetaN: angle between sun vector and facet's normal vector`, scalar (rad)
%[text] `sunB:` sun vector from satellite to Sun expressed with body-fixed frame, 1x3 vector (Nx3 allowed)
%[text] `d`: distance between satellite and Sun, m
%[text] `const`: constants for orbital propagation
%[text] `srp`: approximated SRP with Cook–Torrance model **expressed with Sun-fixed frame**
%[text] ## note
%[text] NA
%[text] ## references 					
%[text] Analytic Approximation of High-Fidelity Solar Radiation Pressure.
%[text] ## revisions
%[text] 20200915  y.yoshimura, y.yoshimula@gmail.com
%[text] See also srpApproxCT2, ctM.
function srp = srpApproxCT(sat, thetaN, sunB, d, const)
arguments (Input)
    sat
    thetaN (:,1)
    sunB (:,3) {mustBeNumeric}
    d (:,1) {mustBeNumeric}
    const
end
arguments (Output)
    srp
end
%[text] ### coefficient
dAU = km2AU(d ./ 10^3, const); % AU
S0 = const.S0; % Solar constant, W/m^2
c = const.c; % light speed, m/s
coeff = -S0 / c / dAU^2;

sat.normal = [0.0, sin(thetaN), cos(thetaN)]; % normal vector
%[text] ### perfect mirror-like reflecttion vector
rRef = 2 * (sunB * sat.normal') .* sat.normal - sunB; % perfect specular direction, 1x3 vector
M = ctM(sat, rRef, sunB); % remaining term
lam = 2 / sat.mCT^2;
c = 1;

sunlitFlag = (sat.normal * sunB' > 0); % nFacetx1 matrix, 1: sunlit, 0: shade

%[text] #### quater-sphere, +y方向の積分範囲
nTheta = 6; % theta分割数
nPhi = 6;
thetaW = (pi/2+thetaN) / 2 / nTheta; % theta_width
phiW = pi / nPhi;
[phiBound, thetaBound] = meshgrid(0:phiW:pi, 0:thetaW:(pi/2+thetaN)/2);
[alp_, bet_] = calcCoeff(phiBound, thetaBound, thetaN, sat.normal, lam, c);

A = zeros(1,3);
B = zeros(1,3);
for j = 1:size(phiBound,1)-1
    for k = 1:size(phiBound,2)-1
        [Atmp, Btmp] = analyticSolCT(thetaN, alp_(j,k), bet_(j,k), ...
            [phiBound(j,k), phiBound(j,k+1)], [thetaBound(j,k), thetaBound(j+1,k)]);
        A = A + Atmp;
        B = B + Btmp;
    end
end

%[text] #### partial hemisphere, -y方向の積分範囲
thetaW = (pi/2 - thetaN)/ 2 / nTheta; % theta_width
phiW = pi / nPhi;
[phiBound, thetaBound] = meshgrid(pi:phiW:2*pi, 0:thetaW:(pi/2-thetaN)/2);
[alp_, bet_] = calcCoeff(phiBound, thetaBound, thetaN, sat.normal, lam, c);
for j = 1:size(phiBound,1)-1
    for k = 1:size(phiBound,2)-1
        [Atmp, Btmp] = analyticSolCT(thetaN, alp_(j,k), bet_(j,k), ...
            [phiBound(j,k), phiBound(j,k+1)], [thetaBound(j,k), thetaBound(j+1,k)]);
        A = A + Atmp;
        B = B + Btmp;
    end
end
srp = sunlitFlag .* coeff .* sat.area .* M .* (A + B);

end

%[appendix]{"version":"1.0"}
%---
