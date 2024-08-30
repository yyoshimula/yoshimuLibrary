function [sat, cd, cs] = lcSimple(sat, sunB, obsB, thr)
% ----------------------------------------------------------------------
%   calculate Light curves with Lambertian and Perfect Specular model
%    20200909  y.yoshimura
%    Inputs: sat, satellite shape definition by obj files, N facets
%            sunB, sun vector from satellite to sun, unit vector, Mx3 matrix 
%            obsB, observer vector from spacecraft to observer, unit vetor, Mx3 
%            thr(optional): perfect specular threshold, rad
%   Outputs: sat: sat.fObs is added, NxM matrix
%   related function files:
%   note:
%   cf:
%   revisions;
% function sat = simpleBRDF(sat, sun_b, obs_b, thr)
%   (c) 2022 yasuhiro yoshimura
%----------------------------------------------------------------------
switch nargin %入力引数の数で場合わけ
    case 3
        thr = deg2rad(1); % perfect specular threshold
    case 4
        thr = thr;
    otherwise
        error('invalid number of arguments')
end
sunB = sunB ./ vecnorm(sunB, 2, 2); %一応 normalize
obsB = obsB ./ vecnorm(obsB, 2, 2);

h = sunB + obsB; % bisector vector of sun and observer vectors
h = h ./ vecnorm(h, 2, 2); % Mx3 matrix

NS = sat.normal * sunB'; % NxM
NV = sat.normal * obsB'; % NxM

% sat.normal * h'は，faceの数 N x 時間履歴の数 Mのmatrixになる
cTotal = sat.Cd./pi ...
    + 2.0 .* sat.Cs .* ((sat.normal * h') >= cos(thr)); 
tmp = cTotal .* sat.area .* (sat.normal * sunB') .* (sat.normal * obsB'); 

% if the faces can be seen or not
cd = sat.Cd ./ pi .* (NS > 0) .* (NV > 0);
cs = 2.0 .* sat.Cs .* ((sat.normal * h') >= cos(thr)) .* (NS > 0) .* (NV > 0); 
sat.fObs = tmp .* (NS > 0) .* (NV > 0); % NxM

end
