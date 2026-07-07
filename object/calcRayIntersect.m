%[text] # calculating ray triangle intersection
%[text] 単一の三角facetに対してray-triangle交差判定を行い，j番目のfacetがi番目のfacetに影を作るかを計算
%[text] ## inputs
%[text] `sun`: Sun directional vector, 1x3
%[text] `nJ`: normal vector of the j-th (shadowing) facet, 1x3
%[text] `vertJ`: vertices of the j-th (shadowing) facet, (x, y, z), m, 3x3 matrix (3頂点)
%[text] `vertI`: position of the i-th (shadowed) facet, m, 1x3
%[text] ## outputs
%[text] `flag`: intersection flag, 1: not shadowed, 0: shadowed
%[text] ## note
%[text] calculate shadow if j-th face make the shadow on the i-th face
%[text] $j$番目のfacetが$i$番目のfacetに影を作るかを判定
%[text] ## references
%[text] NA
%[text] ## revisions
%[text] 20241125 added arguments y.yoshimura
%[text] 20200811  y.yoshimura, y.yoshimula@gmail.com
%[text] See also showSC
function flag = calcRayIntersect(sun, nJ, vertJ, vertI) %#codegen
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
