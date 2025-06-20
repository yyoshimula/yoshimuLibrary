%[text] # calculating object facet area
%[text] 読み込んだ衛星形状の表面積を計算
%[text] `sat`: 構造体satに各種データを格納
%[text] `sat.area`: face area, ${\\rm m^2}$, Nx1 vector, 面の面積
%[text] `sat.pos`: center of face, m, Nx3 matrix, 面の中心（平均）座標
%[text] `sat.vertices`: vertex position, (x, y, z), m, Nx3 matrix, 面を構成する点の座標
%[text] `sat.faces`: face indices, Nx3 matrix, 面を構成する座標のindex
%[text] `sat.Ca`: coefficients for absorption, Nx1 vector, 吸収率
%[text] `sat.Cd`: coefficients for diffusion, Nx1 vector, 拡散反射率
%[text] `sat.Cs`: coefficients for specular reflection, Nx1 vector, 鏡面反射率
%[text] `sat.n`: normal vector, Nx3 matrix, 法線方向単位ベクトル, 外向きが正
%[text] `sat.shadowFlag`: self shadow flag, (default)1: not shadowoed, 0: shadowed
%[text] ## note
%[text] NA
%[text] ## references
%[text] NA
%[text] ## revisions
%[text] 20250120 y.yoshimura 3角メッシュと4角メッシュが両方ある場合に対応
%[text] 2021020209  y.yoshimura
%[text] See also showSC, readSC.
function [area, pos] = calcAreaObj(sat)

area = zeros(size(sat.faces, 1), 1);
pos = zeros(size(sat.faces, 1), 3);

for i = 1:size(sat.faces, 1)

    if size(sat.faces(i,:), 2) == 3 | isnan(sat.faces(i,4))
        nPolygon = 3;
    elseif size(sat.faces(i,:), 2) == ~isnan(sat.faces(i,4))
        nPolygon = 4;
    end

    % vector from index 1 to index 2, nx3 matrix
    vA = sat.vertices(sat.faces(i,2),:) - sat.vertices(sat.faces(i,1),:);

    % vector from index 1 to index 3, nx3 matrix
    vB = sat.vertices(sat.faces(i,3),:) - sat.vertices(sat.faces(i,1),:);

    crossA = cross(vA, vB); % cross product between vA and vB

    if nPolygon == 3
        posSum = sat.vertices(sat.faces(i,1),:) + ...
            sat.vertices(sat.faces(i,2),:) + ...
            sat.vertices(sat.faces(i,3),:); % nx3 matrix

        pos(i,:) = posSum ./ 3.0;
        area(i,1) = vecnorm(crossA, 2, 2) ./ 2; % (norm of cross product) ./ 2

    elseif nPolygon == 4

        posSum = sat.vertices(sat.faces(i,1),:) + ...
            sat.vertices(sat.faces(i,2),:) + ...
            sat.vertices(sat.faces(i,3),:) + ...
            sat.vertices(sat.faces(i,4),:); % nx3 matrix

        pos(i,:) = posSum ./ 4.0;
        area(i,1) = vecnorm(crossA, 2, 2); % (norm of cross product)

    else
        error('polygon must be triangular or quad')
    end
end

end

%[appendix]{"version":"1.0"}
%---
