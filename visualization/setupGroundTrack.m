%% setupGroundTrack - 地上軌跡プロットの背景設定
%
% 地球テクスチャを背景に表示し、軸を設定する

function setupGroundTrack()

libPath = '/Users/yyoshimula/Dropbox/MATLAB/yoshimuLibrary';

earthImg = imread(fullfile(libPath, 'visualization', 'naturalEarth.jpg'));
imagesc([-180, 180], [90, -90], earthImg);
set(gca, 'YDir', 'normal');
hold on;

xlabel('経度 [deg]', 'Color', 'w');
ylabel('緯度 [deg]', 'Color', 'w');
xlim([-180, 180]);
ylim([-90, 90]);
daspect([1 1 1]);

ax = gca;
ax.XColor = 'w';
ax.YColor = 'w';

end
