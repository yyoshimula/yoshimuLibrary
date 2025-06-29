function sat = calcNormal_obj(sat)
% ----------------------------------------------------------------------
%   calculate (outward) normal vectors for .obj format
%    20200903  y.yoshimura
%    Inputs: sat
%   Outputs: sat: sat構造体にnormalを付加して出力
%   related function files:
%   note:
%   cf:
%   revisions;
% function sat = calcNormalObj(sat)
%   (c) 2020 yasuhiro yoshimura
%----------------------------------------------------------------------

% vector from index 1 to index 2, nx3 matrix
v1 = sat.vertices(sat.faces(:,2),:) - sat.vertices(sat.faces(:,1),:);

% vector from index 1 to index 3, nx3 matrix
v2 = sat.vertices(sat.faces(:,3),:) - sat.vertices(sat.faces(:,1),:);


% 	/* 外積 cross productを求める */
crossV = cross(v1, v2);

% 	/* 外積v2×v1の長さ|v2×v1|（= length）を求める */
crossV_norm = vecnorm(crossV, 2, 2);

% 	/* 長さ|v2×v1|が0のときは法線ベクトルは求められない */
crossV_norm = crossV_norm + (crossV_norm <= 1e-8) .*1e-8; % 0割を避けるため微小項をたす

% 	/* 外積v2×v1を長さ|v2×v1|で割って法線単位ベクトルnを求める *
n = crossV ./ crossV_norm;

sat.normal = n;

end