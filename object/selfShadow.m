%[text] # calculate self-shadowing
%[text] ## inputs
%[text] `sat`: satellite configuration
%[text] `sun`: Sun directional vector w.r.t. body-fixed frame, 1x3
%[text] ## note
%[text] 全componentのメッシュに対して処理を行う 
%[text] 影の処理関数($j\n$番目のfacetが影を作る方，$i$番目のfacetは影が写る方)
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20241125 と四角facetに対応 y.yoshimura
%[text] 20200811  y.yoshimura, y.yoshimula@gmail.com
%[text] See also selfShadow.
function sat = selfShadow(sat, sun) %#codegen

sun = sun(:); % column vector
sun = sun ./ norm(sun);

% 影判定用フラグ: 0のとき，影になっているとする
sat.sunlitFlag = zeros(length(sat.faces), 1); % initialization
%%
%[text] ## self-shadowing
% 両方のfacetが太陽光を受けていたら計算する
sunlitIndex = (sat.normal * sun) > 0.0;

faces = sat.faces(sunlitIndex,:);
normal = sat.normal(sunlitIndex,:);
pos = sat.pos(sunlitIndex,:);

tmpFlag = ones(length(faces), 1); % initialization
%[text] use "for", if parallel computation toolbox is not installed
parfor i = 1:length(faces)
    tmpFlagI = 1; % 各iごとの一時フラグ
    posI = pos(i,:);
    for j = 1:length(faces)
        if (j ~= i)
            vertJ = [sat.vertices(faces(j,1),:)
                sat.vertices(faces(j,2),:)
                sat.vertices(faces(j,3),:)]; % 3x3, coordinate of j-th facet
            flag = calcSelfShadow(sun', normal(j,:), vertJ, posI);
            tmpFlagI = tmpFlagI * flag; % OR演算っぽく

            % 四角facetでは2回計算
            if ~isnan(faces(1,4))
                vertJ = [sat.vertices(faces(j,3),:)
                    sat.vertices(faces(j,4),:)
                    sat.vertices(faces(j,1),:)]; % 3x3, coordinate of j-th facet
                flag = calcSelfShadow(sun', normal(j,:), vertJ, posI);
                tmpFlagI = tmpFlagI * flag; % OR演算っぽく
            end

        else
            % do nothing when i == j
        end
    end
    tmpFlag(i,1) = tmpFlagI; % 結果を保存
end

sat.sunlitFlag(sunlitIndex) = tmpFlag;

end


%[appendix]{"version":"1.0"}
%---
