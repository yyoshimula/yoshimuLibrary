%[text] # grid arrangement of figures
%[text] figuresを並べて表示するだけの関数
%[text] `nDisp:` indicating display where figures are arranged,  e.g., 1: primary display, 2: secondary display, etc.
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20240725  y.yoshimura, y.yoshimula@gmail.com
%[text] See also .
function gFigs(nDisp)
arguments
    nDisp (1,1) {mustBeNumeric} = 1
end
wMargin = -50;              % width Margin
hMargin = -100;               % height Margin

% javaを使って各ディスプレイのサイズを取得
screens = java.awt.GraphicsEnvironment.getLocalGraphicsEnvironment().getScreenDevices();

% 指定されたディスプレイのサイズと位置を取得
if nDisp > length(screens)
    error('指定されたディスプレイ番号が存在しません。');
end

% diplayサイズと位置を取得
bounds = screens(nDisp).getDefaultConfiguration().getBounds();
screenWidth = bounds.getWidth();
screenHeight = bounds.getHeight();
screenX = bounds.getX();
screenY = bounds.getY();

figHandles = findobj('Type', 'figure');
nFigs = length(figHandles);

% Define the number of rows and columns
numRows = ceil(sqrt(nFigs));
numCols = ceil(nFigs / numRows);

% Calculate figure width and height
figureWidth = screenWidth / numCols;
figureHeight = screenHeight / numRows;


% Create and position each figure
for i = 1:nFigs

    % 行と列を計算
    row = ceil(i / numCols);
    col = mod(i-1, numCols) + 1;

    % 位置を計算
    left = screenX + (col-1) * figureWidth;
    bottom = screenY + screenHeight - row * figureHeight;

    % Figureの位置を設定
    set(figure(i), 'Position', [left, bottom, wMargin+figureWidth, hMargin+figureHeight]);


end
end

%[appendix]{"version":"1.0"}
%---
