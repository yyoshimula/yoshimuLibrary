%[text] # approximating specular term of SRP with Cook–Torrance model (multi-facet)
%[text] Cook–Torranceモデルを用いたSRP近似解析解（複数 facet 版）。facet ごとに寄与を計算して総和する。
%[text] `sat`: satellite configuration read with `readSC`
%[text] `thetaN: angle between sun vector and each facet's normal vector`, nFacet x 1 column vector (rad)
%[text] `sunB:` sun vector from satellite to Sun expressed with body-fixed frame, 1x3 vector
%[text] `d`: distance between satellite and Sun, m
%[text] `const`: constants for orbital propagation
%[text] `srp`: approximated SRP with Cook–Torrance model **expressed with Sun-fixed frame**, summed over facets
%[text] ## note
%[text] 単一 facet 版 srpApproxCT と異なり、sunlitFlag は 1 固定（日陰判定は無効）。
%[text] ## references
%[text] Analytic Approximation of High-Fidelity Solar Radiation Pressure.
%[text] ## revisions
%[text] 20200915  y.yoshimura, y.yoshimula@gmail.com
%[text] See also srpApproxCT, ctM2.
function srp = srpApproxCT2(sat, thetaN, sunB, d, const)
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
nFacet = size(sat.area,1);
%[text] ### coefficient
dAU = km2AU(d ./ 10^3, const); % AU
S0 = const.S0; % Solar constant, W/m^2
c = const.c; % light speed, m/s
coeff = -S0 / c / dAU^2;

sat.normal = [zeros(nFacet,1), sin(thetaN), cos(thetaN)]; % normal vector, nFacet x 3 matrix
%[text] ### perfect mirror-like reflecttion vector
rRef = 2 * (sat.normal * sunB') .* sat.normal - sunB; % perfect specular direction, nFacet x 3 matrix
M = ctM2(sat, rRef, sunB); % remaining term, nFacet x 1
lam = 2 / sat.mCT.^2;
mu = 1;

% sunlitFlag = (sat.normal * sunB' > 0); % nFacetx1 matrix, 1: sunlit, 0: shade
sunlitFlag = 1;

srpTmp = zeros(nFacet, 3);
for i = 1:nFacet
    satN = sat.normal(i,:);
    lamTmp = lam(i);
%[text] #### quater-sphere, +y方向の積分範囲
    nTheta = 6; % theta分割数
    nPhi = 6;
    thetaW = (pi/2+thetaN(i)) / 2 / nTheta; % theta_width
    phiW = pi / nPhi;
    [phiBound, thetaBound] = meshgrid(0:phiW:pi, 0:thetaW:(pi/2+thetaN(i))/2);
    [alp_, bet_] = calcCoeff(phiBound, thetaBound, thetaN(i), satN, lamTmp, mu);
    
    A = zeros(1,3);
    B = zeros(1,3);
    for j = 1:size(phiBound,1)-1
        for k = 1:size(phiBound,2)-1
            [Atmp, Btmp] = analyticSolCT(thetaN(i), alp_(j,k), bet_(j,k), ...
                [phiBound(j,k), phiBound(j,k+1)], [thetaBound(j,k), thetaBound(j+1,k)]);
            A = A + Atmp;
            B = B + Btmp;
        end
    end

%[text] #### partial hemisphere, -y方向の積分範囲
    thetaW = (pi/2 - thetaN(i))/ 2 / nTheta; % theta_width
    phiW = pi / nPhi;
    [phiBound, thetaBound] = meshgrid(pi:phiW:2*pi, 0:thetaW:(pi/2-thetaN(i))/2);
    [alp_, bet_] = calcCoeff(phiBound, thetaBound, thetaN(i), satN, lamTmp, mu);
    for j = 1:size(phiBound,1)-1
        for k = 1:size(phiBound,2)-1
            [Atmp, Btmp] = analyticSolCT(thetaN(i), alp_(j,k), bet_(j,k), ...
                [phiBound(j,k), phiBound(j,k+1)], [thetaBound(j,k), thetaBound(j+1,k)]);
            A = A + Atmp;
            B = B + Btmp;
        end
    end
    srpTmp(i,:) = sunlitFlag .* coeff .* sat.area(i) .* M(i) .* (A + B);
end
srp = sum(srpTmp,1);
end

%[appendix]{"version":"1.0"}
%---
