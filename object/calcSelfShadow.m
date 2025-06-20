%[text] # calculating self-shadowing of spacecraft
%[text] 読み込んだ衛星形状に対してself-shadowingを計算
%[text] ## note
%[text] calculate shadow if j-th face make the shadow on the i-th face 
%[text] $j\n$番目のfacetが$i\n$番目のfacetに影を作るかを判定
%[text] `satName` : satellite object file name, .obj file
%[text] `sat`: 構造体satに各種データを格納
%[text] `sat.vertices`: vertex position, (x, y, z), m, Nx3 matrix, 面を構成する点の座標
%[text] `sat.faces`: face indices, Nx3 matrix, 面を構成する座標のindex
%[text] `sat.area`: face area, ${\\rm m^2}$, Nx1 vector, 面の面積
%[text] `sat.pos`: center of face, m, Nx3 matrix, 面の中心（平均）座標
%[text] `sat.Ca`: coefficients for absorption, Nx1 vector, 吸収率
%[text] `sat.Cd`: coefficients for diffusion, Nx1 vector, 拡散反射率
%[text] `sat.Cs`: coefficients for specular reflection, Nx1 vector, 鏡面反射率
%[text] `sat.n`: normal vector, Nx3 matrix, 法線方向単位ベクトル, 外向きが正
%[text] `sat.shadowFlag`: self shadow flag, (default)1: not shadowoed, 0: shadowed
%[text] ## references
%[text] NA
%[text] ## revisions
%[text] 20241125 added arguments y.yoshimura
%[text] 20200811  y.yoshimura, y.yoshimula@gmail.com
%[text] See also showSC
function flag = calcSelfShadow(sun, nJ, vertJ, vertI) %#codegen
% arguments
%     sun (1,3) {mustBeNumeric}
%     nJ (:,3) {mustBeNumeric}
%     vertJ (:,3) {mustBeNumeric}
%     vertI (:,3) {mustBeNumeric}
% end

sun = sun(:);
vertI = vertI(:);

d = nJ * vertJ(1,:)'; % scalar

% Calculate intersection scalar K
K = (d - nJ * vertI) / (nJ * sun); % scalar

% 交点計算
Q = vertI + K .* sun; % 3x1
v3 = Q - vertI;
inner = v3' * sun;

%	交点がメッシュ内にあるかどうかの計算
if (inner > 0.0)
    v1 = kron(ones(3,1),Q') - vertJ; % 3x3

    %v2 = 3x3 vector 3点を順番につなぐベクトル
    v2 = [vertJ(2,:) - vertJ(1,:)
          vertJ(3,:) - vertJ(2,:)
          vertJ(1,:) - vertJ(3,:)]; % 3x3 matrix

    % 外積v2×v1（= crossV）を求める //
    crossV = cross(v2, v1, 2); % 3x3
    D = [crossV(1,:) * crossV(2,:)'
        crossV(1,:) * crossV(3,:)'];

    
    if(D(1) >= 0.0 && D(2) >= 0.0)
        flag = 0;
    else
        flag = 1;
    end

else
    flag = 1;
end

end

%[appendix]{"version":"1.0"}
%---
