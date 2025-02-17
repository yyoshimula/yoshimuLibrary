%[text] # calculating light curves using Ashikhmin–Shirley model
%[text] `sat`: satellite configuration read with `readSC, N `facets
%[text] `sunB: sun` vector from satellite, unit vector,  Mx3 matrix
%[text] `obsB`, observer vector from satellite to observer, unit vetor, Mx3 matrix
%[text] `sat: sat.fObs` is added, Nx1
%[text] `cd:` diffuse part of light curves, NxM matrix
%[text] `cs: `specular part of light curves, NxM matrix
%[text] ## note
%[text] 
%[text] ## references 
%[text] Ashikhmin, Michael, & Shirley, Peter. “An Anisotropic Phong BRDF Model.” Journal of graphics tools, vol. 5, no. 2, 2000, pp. 25-32.
%[text] ## revisions
%[text] 20200430  y.yoshimura, y.yoshimula@gmail.com
%[text] See also readSC, lcSimple, lcAS, lcCT, orbitConst.
function [sat, cd, cs] = lcAS(sat, sunB, obsB)
sunB = sunB ./ vecnorm(sunB, 2, 2); %一応 normalize, Mx3
obsB = obsB ./ vecnorm(obsB, 2, 2);
%[text] ### bisector
h = sunB + obsB; % Mx3 matrix, bisector vector of sun and observer vectors
h = h ./ vecnorm(h, 2, 2);

[phiH, thetaH, ~] = cart2sph(h(:,1), h(:,2), h(:,3));
thetaH = pi/2 - thetaH; % +z軸から測った角度にする
phiH = phiH'; % row vector, 1xM
thetaH = thetaH'; % row vector, 1xM
%[text] ### pre-calculation
% sat.normal * sun_b'は，faceの数 N x 時間履歴の数 M のmatrixになる
NS = sat.normal * sunB'; % NxM
NV = sat.normal * obsB'; % NxM
VH = dot(obsB,h,2)'; % 1xM

% reflectance
rho = sat.Cd; % Nx1
F0 = sat.F0;
% NxM matrixへ拡張
nu = repmat(sat.nu, 1, size(sunB,1));
nv = repmat(sat.nv, 1, size(sunB,1));
%[text] ## diffuse
%  NxM matrix
cd = 28 / 23 .* rho / pi .* (1 - F0) .* (1 - (1 - NS./2).^5) .* (1 - (1 - NV./2).^5);
%[text] ## specular
F = F0 + (1 - F0) .* (1 - VH).^5; % NxM
M = sqrt((nu + 1) .* (nv + 1)) ./ 8 ./ pi .* F ./ VH ./ max(NS, NV); % NxM
infIndex = isinf(M);
M(infIndex) = 0; % M計算時にmax関数で0割によるinf発生を除去
D = cos(thetaH).^(nu .* cos(phiH).^2 + nv .* sin(phiH).^2);
D = real(D); % 値が小さすぎて複素数になるのを防ぐため
cs = M .* D; % NxM
%[text] ## total
cTotal = cd + cs;
tmp = cTotal .* sat.area .* NS.* NV;  % NxM
% if the faces can be seen or not
cd = cd .* (NS > 0) .* (NV > 0);
cs = cs .* (NS > 0) .* (NV > 0);
sat.fObs = tmp .* (NS > 0) .* (NV > 0); % NxM
end

%[appendix]{"version":"1.0"}
%---
