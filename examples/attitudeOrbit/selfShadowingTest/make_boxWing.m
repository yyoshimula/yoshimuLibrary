% ----------------------------------------------------------------------
%   definitions of size, split, and optical parameters
%    20190620  y.yoshimura
%    Inputs:
%   Outputs:
%   related function files:
%   note:
%   cf:
%   revisions;
%
%   (c) 2019 yasuhiro yoshimura
%----------------------------------------------------------------------
face = struct;
vert = struct;


%% 光学パラメータ[Ca, Cd, Cs, Ct]
MLI = [0.37; 0.255; 0.375; 0.0];
RAD = [0.1; 0.2; 0.7; 0.0];
% SAP = [0.91; 0.08; 0.01; 0.0];
SAP = MLI;

MAX_EQUIP = 3;
EQUIP_NUM = 6 .* ones(3,1);

%% Satellite Body
body_length = [0.5; 0.5; 0.8];
body_split  = [10; 10; 2];

vert = v4cube(vert, 1, body_length, body_split);
face = f4cube3(face, 1, body_split, MLI);

% the number of facet on each face
MAX_F(1,:) = [body_split(2)*body_split(3), body_split(2)*body_split(3), ...
    body_split(3)*body_split(1), body_split(3)*body_split(1), ...
    body_split(2)*body_split(1), body_split(2)*body_split(1)];

%% SAP
sap_length = [0.5; 1.0; 0.05];
sap_split  = [10; 10; 2];

vert = v4cube(vert, 2, sap_length, sap_split);
face = f4cube3(face, 2, sap_split, SAP);
% vert = rotation(vert, 0, 10, 0, 2, 6); % SAP rotation
vert = translation(vert, 0.0, body_length(2)/2+ sap_length(2)/2.0, 0.0, 2, 6);

MAX_F(2,:) = [sap_split(2)*sap_split(3), sap_split(2)*sap_split(3), ...
    sap_split(3)*sap_split(1), sap_split(3)*sap_split(1), ...
    sap_split(2)*sap_split(1), sap_split(2)*sap_split(1)];
%% SAP
vert = v4cube(vert, 3, sap_length, sap_split);
face = f4cube3(face, 3, sap_split, SAP);
% vert = rotation(vert, 0, 10, 0, 3, 6);
vert = translation(vert, 0.0, -body_length(2)/2-sap_length(2)/2.0, 0.0, 3, 6);

MAX_F(3,:) = [sap_split(2)*sap_split(3), sap_split(2)*sap_split(3), ...
    sap_split(3)*sap_split(1), sap_split(3)*sap_split(1), ...
    sap_split(2)*sap_split(1), sap_split(2)*sap_split(1)];

% for m = 1:MAX_EQUIP
%     for k = 1:EQIP_NUM(m)
%         point{m,k}.coord = point{m,k}.coord - ones(length(point{m,k}.coord),3) * diag(com);    % translation for center of mass
%         point{m,k}.coord = (moi * point{m,k}.coord')'; % rotation for moment of inertia
%     end
% end

MAX_F = MAX_F .* 2.0; % for triangular patch

%% area and normal
face = calcArea3(vert,face);
face = calcNormal(vert, face);
