%[text] # calculating Lambertian diffusion and Perfect Speuclarity
%[text] ## inputs
%[text] `sat`: satellite configuration read with `readSC`
%[text] `sunB:` sun vector from satellite to Sun expressed with body-fixed frame, Mx3 matrix
%[text] `obsB:` observer from satellite with body-fixed frame, Mx3 matrix
%[text] thr: threshold for speuclar lobe, rad
%[text] ## outputs
%[text] fObs: total reflectance BRDF, nFacet x M vector
%[text] cd: diffuse part of BRDF, nFacet x M vector
%[text] cs: specular part of BRDF, nFacet x M vector
%[text] D: normal distribution function of BRDF, nFacet x M vector
%[text] ## note
%[text] NA
%[text] ## references
%[text] NA
%[text] ## revisions
%[text] 20251208  y.yoshimura, y.yoshimula@gmail.com
%[text] See also lcAS, lcCT, readSC.
function [sat, cd, cs] = lcSimple(sat, sunB, obsB, thr)
if nargin < 4
    thr = deg2rad(1); % perfect specular threshold
end
sunB = sunB ./ vecnorm(sunB, 2, 2); %一応 normalize
obsB = obsB ./ vecnorm(obsB, 2, 2);

h = sunB + obsB; % bisector vector of sun and observer vectors
h = h ./ vecnorm(h, 2, 2); % Mx3 matrix

% sat.normal * sunB'は，faceの数 N x 時間履歴の数 Mのmatrixになる
NS = sat.normal * sunB'; % nFacet x M matrix
NV = sat.normal * obsB'; % nFacet x M matrix

cTotal = sat.Cd./pi ...
    + 2.0 .* sat.Cs .* ((sat.normal * h') >= cos(thr));
sat.fObs = cTotal .* sat.area .* NS .* NV;

% if the faces can be seen or not
cd = sat.Cd ./ pi .* (NS > 0) .* (NV > 0);
cs = 2.0 .* sat.Cs .* ((sat.normal * h') >= cos(thr)) .* (NS > 0) .* (NV > 0);
sat.fObs = sat.fObs .* (NS > 0) .* (NV > 0); % nFacet x M matrix

end

%[appendix]{"version":"1.0"}
%---
