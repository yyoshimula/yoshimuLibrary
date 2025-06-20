%[text] # setting and saving figure for presentation slide
%[text] プレゼンテーション用にfigureをいい感じに設定し，保存
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 201XXXXX  y.yoshimura, y.yoshimula@gmail.com
%[text] See also fig4Paper.
function fig4Presen
%[text] ## font 
set(gca,'FontSize',24);
set(gca,'FontName','Times');
set(gca,'FontWeight','normal'); % normal/demi/bold
 
%[text] ## line width
set(gca,'LineWidth', 1.5); % 軸
 
%[text] ## color

clorder = ...
    [  0,   0, 255; % 青
       0, 128,   0; % 緑
     255,   0,   0; % 赤
     204,   8, 204; % 紫
     222, 125,   0; % 茶
       0,  51, 153; % 紺
      64,  64,  64];% 灰色  
set(gca,'ColorOrder',clorder/255);
 
%[text] ## grid
set(gca,'XGrid','on');
set(gca,'YGrid','on');
set(gca,'ZGrid','on');

exportgraphics(gcf,'fig1.pdf', 'ContentType', 'image');

end

%[appendix]{"version":"1.0"}
%---
