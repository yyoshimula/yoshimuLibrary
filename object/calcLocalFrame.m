%[text] # calculating object facet's local frame and its quaternion from body-fixed frame to local frame
%[text] 読み込んだ衛星形状のlocal frameを計算
%[text] ## inputs
%[text] `sat`: 構造体satに各種データを格納
%[text] `sat.normal`: normal vector, Nx3 matrix, 法線方向単位ベクトル, 外向きが正
%[text] `sat.vertices`: vertex position, (x, y, z), m, Nx3 matrix, 面を構成する点の座標
%[text] `sat.faces`: face indices, Nx3 matrix, 面を構成する座標のindex
%[text] ## outputs
%[text] `uu`: x-axis of local frame w.r.t. body-fixed frame, Nx3 matrix
%[text] `uv`: y-axis of local frame w.r.t. body-fixed frame, Nx3 matrix
%[text] `qlb`: quaternion from body-fixed frame to local frame, Nx4 matrix
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
