%[text] # setting and saving figure for conference or journal manuscript
%[text] 論文用にfigureをいい感じに設定し，保存
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20240823  y.yoshimura, y.yoshimula@gmail.com, major update
%[text] See also fig4Presen.
function fig4Paper(fig, nFig, options)
arguments
    fig = gcf
    nFig = 1;
    options.contentType = 'vector'
    options.asis = 0;
end

% colorOrder = ...
%     [216 82 25 % Red
%     118 172 48 % Green
%     0 114 190; % Blue
%     204,   8, 204; % Purple
%     222, 125,   0; % Brown
%     64,  64,  64]; % Dark Gray


% 現在の図サイズを保持（リサイズしない）
fig.PaperPositionMode = 'auto';

%[text] #### find axes & set font
% font setting
adjustFont(fig)

% optimize line width
optimizeFig(fig, options.asis)

if options.asis == 0
    [nRows,nCols] = detectLayout(fig);         % 行列数を自動取得
    [Wcm,Hcm]     = decideFigSize(nRows,nCols);
    set(fig,'Units','centimeters');
    fig.Position(3:4) = [Wcm Hcm];             % 幅・高さのみ更新
else % asis == 1のときはfigure positionを変えない
    % no change of fig.Position
end

% LaTeXラベルを適用
applyLatexLabels(fig);
set(fig,'DefaultTextInterpreter','latex')  % 以降は LaTeX が既定
% set(fig,'DefaultLabelInterpreter','latex')  % 以降は LaTeX が既定
set(fig,'DefaultAxesTickLabelInterpreter','latex')  % 目盛りも統一
set(fig,'DefaultLegendInterpreter','latex')         % 凡例も統一

%[text] #### export
% file name
timeStamp = datetime('now', 'Format', 'yyyyMdd-HHmmss');
fName = strcat('fig', num2str(nFig), '_', string(timeStamp), '.pdf');

exportgraphics(gcf, fName, 'ContentType', options.contentType, 'BackgroundColor', 'none', 'Resolution', 600);

end

%[text] ### フォント設定
function adjustFont(fig)
% AIAA推奨フォント設定（デフォルト図サイズ用）
% Times New Roman, 論文用に適切なサイズ

% デフォルト図サイズに適したフォントサイズ
targetFontSize = 12; % デフォルト図サイズ用

% 軸の設定
axs = findall(fig, 'Type', 'Axes');
for i = 1:length(axs)
    ax = axs(i);

    % フォント設定（フォールバック付き）
    try
        set(ax, 'FontName', 'Times New Roman', 'FontSize', targetFontSize);
    catch
        set(ax, 'FontName', 'Times', 'FontSize', targetFontSize);
    end

    % ラベルのフォントサイズ調整
    if ~isempty(get(ax, 'Title'))
        set(get(ax, 'Title'), 'FontSize', targetFontSize + 2, 'FontWeight', 'bold');
    end
    if ~isempty(get(ax, 'XLabel'))
        set(get(ax, 'XLabel'), 'FontSize', targetFontSize);
    end
    if ~isempty(get(ax, 'YLabel'))
        set(get(ax, 'YLabel'), 'FontSize', targetFontSize);
    end
    if ~isempty(get(ax, 'ZLabel'))
        set(get(ax, 'ZLabel'), 'FontSize', targetFontSize);
    end
end

% テキストオブジェクト
texts = findall(fig, 'Type', 'Text');
for i = 1:length(texts)
    set(texts(i), 'FontSize', targetFontSize, 'FontName', 'Times New Roman');
end

% 凡例
lgs = findall(fig, 'Type', 'Legend');
for i = 1:length(lgs)
    set(lgs(i), 'FontSize', targetFontSize - 1, 'FontName', 'Times New Roman');
end
end

%[text] ### 図の外観最適化
function optimizeFig(fig, asis)
% 図の外観最適化

axs = findall(fig, 'Type', 'Axes');
for i = 1:length(axs)
    ax = axs(i);

    % 線の太さ調整（印刷時の視認性向上）
    lines = findall(ax, 'Type', 'Line');
    for j = 1:length(lines)
        if get(lines(j), 'LineWidth') < 1.0
            set(lines(j), 'LineWidth', 1.0);
        end
    end

    % マーカーサイズ調整
    for j = 1:length(lines)
        markerSize = get(lines(j), 'MarkerSize');
        if markerSize > 0 && markerSize < 6
            set(lines(j), 'MarkerSize', 6);
        end
    end

    % グリッド設定
    set(ax, 'GridLineStyle', '-', 'GridAlpha', 0.3);

    % 軸の設定
    set(ax, 'Box', 'on', 'LineWidth', 0.8);

    % カラーマップ最適化（グレースケール印刷対応）
    if ~isempty(findall(ax, 'Type', 'Surface')) || ~isempty(findall(ax, 'Type', 'Image'))
        colormap(ax, 'gray'); % または適切なカラーマップ
    end

    if asis == 0
        % 凡例の位置最適化
        lgs = findall(fig, 'Type', 'Legend');
        for i = 1:length(lgs)
            set(lgs(i), 'Box', 'on', 'Location', 'best');
        end
    else
        % do nothing
    end
end

end

%[text] ### LaTeXインタープリターを適用
function applyLatexLabels(fig)

axs = findall(fig, 'Type', 'Axes');
for i = 1:length(axs)
    ax = axs(i);

    % 軸ラベルとタイトル
    xlabel_h = get(ax, 'XLabel');
    ylabel_h = get(ax, 'YLabel');
    zlabel_h = get(ax, 'ZLabel');
    title_h = get(ax, 'Title');

    if ~isempty(xlabel_h)
        set(xlabel_h, 'Interpreter', 'latex');
    end
    if ~isempty(ylabel_h)
        set(ylabel_h, 'Interpreter', 'latex');
    end

    if ~isempty(zlabel_h)
        set(zlabel_h, 'Interpreter', 'latex');
    end

    if ~isempty(title_h)
        set(title_h, 'Interpreter', 'latex');
    end
end

% 凡例
lgs = findall(fig, 'Type', 'Legend');
for i = 1:length(lgs)
    set(lgs(i), 'Interpreter', 'latex');
end
end

%[text] ### subplotやtiledlayoutの行数・列数を確認

function [nRows,nCols] = detectLayout(fig)
tl = findall(fig,'Type','TiledChartLayout');
if ~isempty(tl)
    gs         = tl(1).GridSize;
    nRows = gs(1); nCols = gs(2);
else
    ax   = findall(fig,'Type','Axes','-not','Tag','legend');
    pos  = cat(1,ax.Position);
    % 行をグループ化
    yCenters = round(pos(:,2)*100);       % ピクセル誤差回避
    rows     = findgroups(yCenters);
    nRows    = numel(unique(rows));
    % 列をグループ化
    xCenters = round(pos(:,1)*100);
    cols     = findgroups(xCenters);
    nCols    = numel(unique(cols));
end
end

function [Wcm,Hcm] = decideFigSize(nRows,nCols)
% 幅
if nCols==1
    Wcm = 8.3;                        % 1-column
elseif nCols==2 && nRows==1
    Wcm = 17.0;                       % 横並び2枚 → 2-column
else
    Wcm = min(8.3*nCols*0.9,17.0);    % 3枚以上でも上限17 cm
end

% 高さ
baseH = 8.3*0.68;                     % 1段の標準高さ
Hcm   = baseH*nRows*0.95;             % 行数に比例．余白を少し圧縮
end

%[appendix]{"version":"1.0"}
%---
