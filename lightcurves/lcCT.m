function [sat, cd, cs] = lcCT(sat, sunB, obsB, NDF) 
% ----------------------------------------------------------------------
%   calculate Light curves with Cook–Torrance model
%    20220317  y.yoshimura
%    latest update: 20220317
%    Inputs: sat, satellite shape definition by obj files, N facets
%            sunB, sun vector from satellite to sun, unit vector, Mx3 matrix 
%            obsB, observer vector from spacecraft to observer, unit vetor, Mx3             
%   Outputs: cd: diffuse, Nx1
%            cs: specualr, Nx1
%            sat: sat.fObs is added, Nx1
%   related function files:
%   note:
%   cf:
%   revisions;
%   function [sat, cd, cs] = lcCT(sat, sunB, obsB) 
%   (c) 2022 yasuhiro yoshimura
%----------------------------------------------------------------------
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

sunB = sunB ./ vecnorm(sunB, 2, 2); %一応 normalize, Mx3
obsB = obsB ./ vecnorm(obsB, 2, 2);

h = sunB + obsB; % Mx3 matrix, bisector vector of sun and observer vectors
h = h ./ vecnorm(h, 2, 2);
[phi_h, theta_h, ~] = cart2sph(h(:,1), h(:,2), h(:,3));
theta_h = pi/2 - theta_h; % +z軸から測った角度にする
phi_h = phi_h'; % row vector, 1xM
theta_h = theta_h'; % row vector, 1xM

% reflectance
rho = sat.Cd; % Nx1
F0 = sat.Cs;
m = sat.m;

% sat.normal * sun_b'は，faceの数 N x 時間履歴の数 M のmatrixになる
NS = sat.normal * sunB'; % NxM
NH = sat.normal * h'; % NxM
NV = sat.normal * obsB'; % NxM
VH = dot(obsB,h,2)'; % 1xM

%% diffuse
cd = rho ./ pi;

%% specular
nest = (1 + sqrt(F0)) ./ (1 - sqrt(F0));
g = sqrt(nest.^2 + VH.^2 - 1); % NxM

if strcmp(NDF, 'Beckmann')% Beckmann distribution
    D = exp(-(tan(theta_h)./m).^2); % NxM
    D = D ./ pi ./ m.^2 ./ cos(theta_h).^4;
else
    D = exp(-(theta_h./m).^2); % Gaussian distribution
end

temp1 =  2 * NH .* NV ./ VH; % NxM
temp2 =  2 * NH .* NS ./ VH;

G = min(1, temp1);
G = min(G, temp2);

temp1 = (g - VH).^2 ./ 2 ./ (g + VH).^2;
temp2 = (1 + (VH .* (g + VH) - 1).^2 ./ (VH * (g - VH) + 1).^2);
F = temp1 * temp2; % NxM

cs = D .* G .* F ./ (NS) ./ (NV) ./ 4;

tmp = (cd + cs) .* sat.area .* NS.* NV;  % NxM

% if the faces can be seen or not
cd = cd .* (NS > 0) .* (NV > 0);
cs = cs .* (NS > 0) .* (NV > 0);
sat.fObs = tmp .* (NS > 0) .* (NV > 0); % NxM
end