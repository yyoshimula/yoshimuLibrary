%[text] # setting and saving figure for conference or journal manuscript
%[text] 論文用にfigureをいい感じに設定し，保存
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20240823  y.yoshimura, y.yoshimula@gmail.com, major update
%[text] See also fig4Presen.
function fig4Paper(contentType)
arguments
    contentType = 'vector'
end

clorder = ...
    [216 82 25 % Red
    118 172 48 % Green
    0 114 190; % Blue
    204,   8, 204; % Purple
    222, 125,   0; % Brown
    64,  64,  64]; % Dark Gray
%[text] ## find fig & its axes
fig = gcf;
fig.PaperType = 'a4';
fig.PaperUnits      = 'centimeters';
fig.Units           = 'centimeters';

% Define the position and size of the paper
fig.PaperPosition   = [2,2, 4, 4];
% fig.Position        = [2,2, 6, 6]; % Uncomment to set figure window size

% Set the figure background color to white
fig.Color = 'w';
% Prevent the figure from inverting colors when printing
fig.InvertHardcopy  = 'off';
% Automatically retrieve all axes objects within the figure
ax = findobj(fig, 'Type', 'Axes');

%[text] ## set parameters
% 取得した軸オブジェクトに対して操作を行う
for i = 1:length(ax)

%[text] ### font 
    ax(i).FontSize = 20;
    ax(i).FontName = 'Times New Roman';
    ax(i).FontWeight = 'normal';

%[text] ### line width
    ax(i).LineWidth = 1.0; % 軸

%[text] ### color
    ax(i).ColorOrder = clorder/255;

%[text] ### grid
    ax(i).XGrid = 'off';
    ax(i).YGrid = 'off';
    ax(i).ZGrid = 'off';

end

%[text] ## export
% file name
timeStamp = datetime('now', 'Format', 'yyyyMdd-HHmmss');
fName = strcat(string(timeStamp), '.pdf');

exportgraphics(gcf, fName, 'ContentType', contentType);

end

%[appendix]{"version":"1.0"}
%---
