% calculate self-shadowing
% inputs
% note
% 全componentのメッシュに対して処理を行う
% 影の処理関数(j番目のfacetが影を作る方，i番目のfacetは影が写る方)
% references
% NA
% revisions
% 20241125 三角facetと四角facetの両方に対応
% 20200811  y.yoshimura, y.yoshimula@gmail.com
% See also q2zyx.

function sat = selfShadow_(sat, sun) %#codegen
sun = sun(:); % column vector
sun = sun ./ norm(sun);

% 影判定用フラグ: 0のとき，影になっているとする
sat.shadowFlag = zeros(length(sat.faces), 1); % initialization

% 両方のfacetが太陽光を受けていたら計算する
sunlitIndex = (sat.normal * sun) > 0.0;

faces = sat.faces(sunlitIndex,:);
normal = sat.normal(sunlitIndex,:);
pos = sat.pos(sunlitIndex,:);

tmpFlag = ones(length(faces), 1); % initialization

% use "for", if parallel toolbox is not installed
parfor i = 1:length(faces)
    tmpFlagI = 1; % 各iごとのtemporaryフラグ
    posI = pos(i,:);
    for j = 1:length(faces)
        if (j ~= i)

            vertJ = [sat.vertices(faces(j,1),:)
                sat.vertices(faces(j,2),:)
                sat.vertices(faces(j,3),:)]; % 3x3, coordinate of j-th facet
            flag = calcSelfShadow_(sun, normal(j,:), vertJ, posI);
            tmpFlagI = tmpFlagI * flag; % OR演算っぽく

            % 四角facetでは2回計算 
            if size(faces,2) == 4
                vertJ = [sat.vertices(faces(j,3),:)
                    sat.vertices(faces(j,4),:)
                    sat.vertices(faces(j,1),:)] % 3x3, coordinate of j-th facet
                flag = calcSelfShadow_(sun, normal(j,:), vertJ, posI);
                tmpFlagI = tmpFlagI * flag; % OR演算っぽく
            end

        else
            % do nothing when i == j
        end
    end
    tmpFlag(i,1) = tmpFlagI; % 結果を保存
end

sat.shadowFlag(sunlitIndex) = tmpFlag;

end
