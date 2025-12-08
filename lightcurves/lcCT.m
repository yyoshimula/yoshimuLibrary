%[text] # calculating Cook–Torrance model
%[text] ## inputs
%[text] `sat`: satellite configuration read with `readSC`
%[text] `sunB:` sun vector from satellite to Sun expressed with body-fixed frame, Mx3 matrix
%[text] `obsB:` observer from satellite with body-fixed frame, Mx3 matrix
%[text] NDF: NDF definition (option), 'Beckamnn' (default) or 'Gauss'
%[text] ## outputs
%[text] fObs: total reflectance BRDF, nFacet x M vector
%[text] cd: diffuse part of BRDF, nFacet x M vector
%[text] cs: specular part of BRDF, nFacet x M vector
%[text] D: normal distribution function of BRDF, nFacet x M vector
%[text] ## note
%[text] Cook–Torrance model is written as
%[text] $c\_{d} = \\frac{\\rho\_{d}}{\\pi}$
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
%[text] 20251208  y.yoshimura, y.yoshimula@gmail.com
%[text] See also lcAS, readSC.
function [sat, cd, cs, D] = lcCT(sat, sunB, obsB, NDF)
switch nargin %入力引数の数で場合わけ
    case 3
        NDF = 'Beckmann'; % NDFのdistribution, default
    case 4
        if strcmp(NDF, 'Gauss')
            NDF = 'Gauss'; %最後の引数にGaussという文字列を入れた場合Gaussian NDF
        else
            NDF = 'Beckmann';
        end
    otherwise
        NDF = 'Beckmann';
end

M = size(sunB, 1);
sunB = sunB ./ vecnorm(sunB, 2, 2); %一応 normalize, Mx3
obsB = obsB ./ vecnorm(obsB, 2, 2);
h = sunB + obsB; % Mx3 matrix, bisector vector of sun and observer vectors
h = h ./ vecnorm(h, 2, 2);

% sat.normal * sunB'は，nFacet x M のmatrixになる
NS = sat.normal * sunB'; % nFacet x M
NH = sat.normal * h'; % nFacet x M
NV = sat.normal * obsB'; % nFacet x M
VH = dot(obsB,h,2)'; % 1xM

thetaH = acos(NH); % nFacet x M
%%
%[text] ## diffuse
cd = sat.Cd ./ pi;
cd = repmat(cd, 1, M);
%%
%[text] ## specular
nest = (1 + sqrt(sat.F0)) ./ (1 - sqrt(sat.F0));
g = sqrt(nest.^2 + VH.^2 - 1); % nFacet x M

if strcmp(NDF, 'Beckmann')% Beckmann distribution
    D = exp(-(tan(thetaH)./sat.mCT).^2); % nFacet x M
    D = D ./ pi ./ sat.mCT.^2 ./ cos(thetaH).^4;
elseif strcmp(NDF, 'Gauss')
    D = exp(-(thetaH./sat.mCT).^2); % Gaussian distribution
else
    error('set the proper NDF option')
end

temp1 =  2 * NH .* NV ./ VH; % nFacet x M
temp2 =  2 * NH .* NS ./ VH;

G = min(1, temp1);
G = min(G, temp2);

temp1 = (g - VH).^2 ./ 2 ./ (g + VH).^2;
temp2 = (1 + (VH .* (g + VH) - 1).^2 ./ (VH * (g - VH) + 1).^2);
F = temp1 * temp2; % nFacet x M

cs = D .* G .* F ./ NS ./ NV ./ 4;

tmp = (cd + cs) .* sat.area .* NS.* NV;  % nFacet x M

% if the faces can be seen or not
cd = cd .* (NS > 0) .* (NV > 0);
cs = cs .* (NS > 0) .* (NV > 0);
sat.fObs = tmp .* (NS > 0) .* (NV > 0); % nFacet x M
end

%[appendix]{"version":"1.0"}
%---
