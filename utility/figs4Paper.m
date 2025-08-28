%[text] # setting and saving figure for conference or journal manuscript
%[text] 論文用にfigureをいい感じに設定し，保存
%[text] 表示されているfig全てを保存
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20240823  y.yoshimura, y.yoshimula@gmail.com, major update
%[text] See also fig4Presen.
function figs4Paper(asis, contentType)
if nargin < 1
    asis = 0;
elseif nargin < 2
    contentType = 'vector';
end

% 図のUserDataプロパティを使って保存状態を記録

%% ❶ すべての可視 Figure を取得
figs = findall(groot,'Type','figure','Visible','on');

%% ❂ Live Script インライン Figure を除外
isInlineLE = startsWith({figs.Tag},'TMWLiveEditor');  % ←最重要フィルタ
figs       = figs(~isInlineLE);

%% ❸ Figure 番号で昇順ソート
[~,ord] = sort([figs.Number]);
figs    = figs(ord);

%% ❹ ループ処理（Live Script対応の重複防止）
for k = 1:numel(figs)
    currentFig = figs(k);
    
    if ~isvalid(currentFig)
        continue;
    end
    
    % UserDataを使って保存済みかチェック
    userData = get(currentFig, 'UserData');
    if ~isstruct(userData) || ~isfield(userData, 'fig4Papers_saved')
        % 未保存の場合
        fig4Paper(asis, currentFig, k, contentType);  % 保存
        
        % 保存済みマークを設定
        if ~isstruct(userData)
            userData = struct();
        end
        userData.fig4Papers_saved = true;
        userData.fig4Papers_timestamp = datetime('now');
        set(currentFig, 'UserData', userData);
        
        fprintf('Saved: Figure %d\n', currentFig.Number);
    else
        fprintf('Skipped (already saved): Figure %d\n', currentFig.Number);
    end
end

lastSaveTime = datetime('now');
end


%[appendix]{"version":"1.0"}
%---
