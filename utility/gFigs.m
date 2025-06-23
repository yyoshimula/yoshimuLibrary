%[text] # grid arrangement of figures
%[text] figuresを並べて表示するだけの関数
%[text] `nDisp:` indicating display where figures are arranged,  e.g., 1: primary display, 2: secondary display, etc.
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20250623  y.yoshimura, y.yoshimula@gmail.com
%[text] See also .
function gFigs(nDisp)
% gFigs  ― 既存の Figure を選んだディスプレイ上にグリッド配置する
%
%  gFigs()          : プライマリディスプレイ（1）へ配置
%  gFigs(nDisp)     : nDisp で指定したディスプレイへ配置（1,2,…）
%
%  2025-07-25  y.yoshimura
%
arguments
    nDisp (1,1) {mustBeInteger,mustBePositive} = 1
end

%%--- ディスプレイ情報
monPos = get(0,'MonitorPositions');   % [left bottom width height]（pixel）
if nDisp > size(monPos,1)
    error('指定されたディスプレイ番号が存在しません．');
end
scr = monPos(nDisp,:);                % 対象ディスプレイ

%%--- Figure ハンドル取得
figs = findall(0,'Type','figure','Visible','on'); % Visible な Figure
nF   = numel(figs);
if nF == 0,  warning('表示中の Figure がありません．');  return;  end

%%--- グリッドサイズ計算
nRow = ceil(sqrt(nF));
nCol = ceil(nF/nRow);

gap   = 50;                                   % マージン（ピクセル）
fW    = (scr(3) - gap*(nCol+1)) / nCol;       % 幅
fH    = (scr(4) - gap*(nRow+1)) / nRow;       % 高さ

%%--- 各 Figure を移動
for k = 1:nF
    r = ceil(k/nCol);              % 行（1-based）
    c = k - (r-1)*nCol;            % 列

    left   = scr(1) + gap + (c-1)*(fW+gap);
    bottom = scr(2) + scr(4) - r*(fH+gap);

    set(figs(k), 'Units','pixels', 'Position',[left bottom fW fH]);
end
end

%[appendix]{"version":"1.0"}
%---
