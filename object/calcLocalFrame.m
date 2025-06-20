%[text] # calculating object facet's local frame and its quaternion from body-fixed frame to local frame
%[text] 読み込んだ衛星形状のlocal frameを計算
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
%[text] 20230828  y.yoshimura
%[text] See also showSC, readSC, calcAreaObj.
function [uu, uv ,qlb] = calcLocalFrame(sat)

% size of satellite's face
n = size(sat.faces, 1);
qlb = zeros(n,4);

% vector from index 1 to index 2, nx3 matrix
v1 = sat.vertices(sat.faces(:,2),:) - sat.vertices(sat.faces(:,1),:);

uu = v1 ./ vecnorm(v1 , 2, 2); % x-axis of local frame expressed with body-fixed frame


% directional cosine matrix using triad method
dcm = triad(uu, sat.normal, repmat([1 0 0], n, 1), repmat([0 0 1], n, 1));

% quaternion
for i = 1:n
    qlb(i,:) = dcm2q(4, dcm(:,:,i));
end

uv = cross(sat.normal, uu);


end


%[appendix]{"version":"1.0"}
%---
