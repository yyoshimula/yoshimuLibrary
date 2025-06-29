function sat = calcShadow_obj(sat, sun)
% ----------------------------------------------------------------------
%   calculate self-hadowin
%    20200903  y.yoshimura
%    Inputs:
%   Outputs:
%   related function files:
%   note:
%   cf:
%   revisions;
%
%   (c) 2020 yasuhiro yoshimura
%----------------------------------------------------------------------
% 全componentのメッシュに対して処理を行う
% //■影の処理関数(j番目のfaceが影を作るほう，i番目のfaceが影が写るほうのメッシュ)

sun = sun(:); % column vector

% 影判定用フラグ 0のとき，影になっているとする
sat.shadowFlag = ones(length(sat.faces), 1); % initialization

%% self-shadowing
for i = 1:length(sat.faces)
    for j = 1:length(sat.faces)
        if (j ~= i)
            vert_j = [sat.vertices(sat.faces(j,1),:)
                sat.vertices(sat.faces(j,2),:)
                sat.vertices(sat.faces(j,3),:)]; % 3x3, coordinate of j-th face      
            flag = selfShadow_obj(sun, sat.normal(j,:), vert_j, sat.normal(i,:), sat.pos(i,:));
            sat.shadowFlag(i,1) = sat.shadowFlag(i,1) * flag; % OR演算っぽく
        else
            % do nothing when i == j
        end
    end
end

end