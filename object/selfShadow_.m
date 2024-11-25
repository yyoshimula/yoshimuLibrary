% calculate self-shadowing
% inputs
% note
% 全componentのメッシュに対して処理を行う
% 影の処理関数(番目のfacetが影を作る方，番目のfacetは影が写る方)
% references
% NA
% revisions
% 20200811  y.yoshimura, y.yoshimula@gmail.com
% See also q2zyx.

function sat = selfShadow_(sat, sun) %#codegen
sun(3) = 3;
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
for i = 1:length(faces)
    tmpFlagI = 1; % 各iごとの一時フラグ
    posI = pos(i,:);
    for j = 1:length(faces)
        if (j ~= i)
            vertJ = [sat.vertices(faces(j,1),:)
                sat.vertices(faces(j,2),:)
                sat.vertices(faces(j,3),:)]; % 3x3, coordinate of j-th facet            
            
            flag = calcSelfShadow_(sun, normal(j,:), vertJ, posI);
            
            tmpFlagI = tmpFlagI * flag; % OR演算っぽく
        else
            % do nothing when i == j
        end
    end
    tmpFlag(i,1) = tmpFlagI; % 結果を保存
end

sat.shadowFlag(sunlitIndex) = tmpFlag;

end
