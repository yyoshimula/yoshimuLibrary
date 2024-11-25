function flag = calcSelfShadow_(sun, nJ, vertJ, vertI)
% calculate shadow if j-th face make the shadow on the i-th face
% j番目のfacetがi番目のfacetに影を作るかを判定
% 20241125 三角facetと四角facetの両方に対応

sun = sun(:);
vertI = vertI(:); % 3x1 vector

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