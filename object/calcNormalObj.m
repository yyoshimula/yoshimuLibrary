%[text] # calculating object facet normal vectors
%[text] 読み込んだ衛星形状の法線ベクトルを計算
%[text] `sat`: 構造体satに各種データを格納
%[text] `sat.n`: normal vector, Nx3 matrix, 法線方向単位ベクトル, 外向きが正
%[text] `sat.area`: face area, ${\\rm m^2}$, Nx1 vector, 面の面積
%[text] `sat.pos`: center of face, m, Nx3 matrix, 面の中心（平均）座標
%[text] `sat.vertices`: vertex position, (x, y, z), m, Nx3 matrix, 面を構成する点の座標
%[text] `sat.faces`: face indices, Nx3 matrix, 面を構成する座標のindex
%[text] `sat.Ca`: coefficients for absorption, Nx1 vector, 吸収率
%[text] `sat.Cd`: coefficients for diffusion, Nx1 vector, 拡散反射率
%[text] `sat.Cs`: coefficients for specular reflection, Nx1 vector, 鏡面反射率
%[text] `sat.shadowFlag`: self shadow flag, (default)1: not shadowoed, 0: shadowed
%[text] ## note
%[text] NA
%[text] ## references
%[text] NA
%[text] ## revisions
%[text] 2021020209  y.yoshimura
%[text] See also showSC, readSC, calcAreaObj.
function n = calcNormalObj(sat)

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

end

%[appendix]{"version":"1.0"}
%---
