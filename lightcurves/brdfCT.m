%[text] # bidrectional reflectance distribution function (BRDF) of Cook–Torrance model
%[text] Cook–Ashikhmin–Shirleyモデルのdiffuse, specular, NDF計算
%[text] ## inputs
%[text] `sat`: satellite configuration read with `readSC`
%[text]  `sat.normal, sat.Cd, sat.mCT` are necessary
%[text] `sunB:` sun vector from satellite to Sun expressed with body-fixed frame, Mx3 matrix
%[text] `vB: reference` expressed with body-fixed frame, Mx3 matrix
%[text] NDF: NDF definition, 'Beckamnn' (default) or 'Gauss'
%[text] ## outputs
%[text] cd: diffuse part of BRDF, nFacet x 1 vector
%[text] cs: specular part of BRDF, nFacet x 1 vector
%[text] D: normal distribution function of BRDF, nFacet x 1 vector
%[text] ## note
%[text] Cook–Torrance model is written as
%[text] $c\_{d} = \\frac{\\rho}{\\pi}$
%[text] $c\_{s} =  \\frac{DGF}{4\\left({\\bf n}^{T}\\bf{s}\\right)\\left(\\bf{n}^{T}\\bf{v}\\right)}$
%[text] where
%[text] Beckmann distribution: $D = \\frac{1}{\\pi m^{2} \\cos^{4}{\\theta\_{h}}}e^{-\\left(\\frac{\\tan{\\theta\_{h}}}{m}\\right)^{2}}$ 
%[text] or Gaussian model: $D=ce^{-(\\alpha/m)^2$
%[text] $G={\\rm min}\\left\\{1,\\frac{2\\left(\\bf{n}^{T}\\bf{h}\\right)\\left(\\bf{n}^{T}\\bf{v}\\right)}{\\bf{v}^{T}\\bf{h}},\\frac{2\\left(\\bf{n}^{T}\\bf{h}\\right)\\left(\\bf{n}^{T}\\bf{s}\\right)}{\\bf{v}^{T}\\bf{h}}\\right\\}$
%[text] $F = \\frac{\\left(g-\\bf{v}^{T}\\bf{h}\\right)^{2}}{2\\left(g+\\bf{v}^{T}\\bf{h}\\right)^{2}}\\left\\{1+\\frac{\\left\[\\bf{v}^{T}\\bf{h}\\left(g+\\bf{v}^T\\bf{h}\\right)-1\\right\]^{2}}{\\left\[\\bf{v}^{T}\\bf{h}\\left(g-\\bf{v}^T{\\bf h}\\right)+1\\right\]^{2}}\\right\\}$
%[text] $g^2 = n\_{\\rm ref}^2 + ({\\bf n^T h})^2 - 1$
%[text] $n\_{\\rm ref} = \\frac{1+\\sqrt{F\_0}}{1-\\sqrt{F\_0}}$
%[text] ## references 
%[text] Cook, R. L., & Torrance, K. E. (1982). A reflectance model for computer graphics. ACM Transactions on Graphics (TOG), 1, 7-24. 
%[text] ## revisions
%[text] 20221208  y.yoshimura, y.yoshimula@gmail.com
%[text] See also brdfAS, readSC.
function [cd, cs, D] = brdfCT(sat, sunB, vB, NDF) 
% arguments
%     sat
%     sunB (:,3) {mustBeNumeric}
%     vB (:,3) {mustBeNumeric}
%     NDF (1,:) char {mustBeMember(NDF,{'Beckmann', 'Gauss'})} = 'Beckmann'
% end

%[text] ### irradiance and radiance
M = size(sunB, 1);

sunB = sunB ./ vecnorm(sunB, 2, 2); %一応 normalize, Mx3
vB = vB ./ vecnorm(vB, 2, 2);

h = sunB + vB; % Mx3 matrix, bisector vector of sun and observer vectors
h = h ./ vecnorm(h, 2, 2);

[~, thetaH, ~] = cart2sph(h(:,1), h(:,2), h(:,3));
thetaH = pi/2 - thetaH; % +z軸から測った角度にする
% phiH = phiH'; % row vector, 1xM
% thetaH = thetaH'; % row vector, 1xM

%[text] ### dot products
% sat.normal * sunB'は，faceの数 nFacet x 時間履歴の数 M のmatrixになる
NS = sat.normal * sunB'; % nFacetxM
NV = sat.normal * vB'; % nFacetxM
NH = sat.normal * h';
VH = dot(vB,h,2)'; % 1xM

%[text] ## diffuse, nFacet x M
cd = sat.Cd / pi;
cd = repmat(cd, 1, M);
%[text] ## specualr, nFacet x M
nRef = (1 + sqrt(sat.F0)) / (1 - sqrt(sat.F0));
g = sqrt(nRef.^2 + VH.^2 - 1);

if strcmp(NDF, 'Beckmann')
    D = exp(-(tan(thetaH)./sat.mCT).^2); % Beckmann distribution
    D = D ./ pi ./ sat.mCT.^2 ./ cos(thetaH).^4;
elseif strcmp(NDF, 'Gauss')
    D = exp(-(thetaH./sat.mCT).^2); % Gaussian distribution
else
    error('set the proper NDF option')
end

tmp1 =  2 * NH .* NV ./ VH;
tmp2 =  2 * NH .* NS ./ VH;

G = min(1, tmp1);
G = min(G, tmp2);

tmp1 = (g - VH).^2 / 2 / (g + VH).^2;
tmp2 = (1 + (VH .* (g + VH) - 1).^2 ./ (VH .* (g - VH) + 1).^2);
F = tmp1 * tmp2;

cs = D .* G .* F ./ NS ./ NV / 4;

end


%[appendix]{"version":"1.0"}
%---
