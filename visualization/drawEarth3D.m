%% drawEarth3D - 3D地球描画（テクスチャ+影オーバーレイ）
%
% 入力:
%   GMST0: 初期GMST [rad]
%   const: 軌道定数
%   sunDir: 太陽方向単位ベクトル (1x3), 空の場合は影なし
%
% 出力:
%   hEarthTransform: 地球回転用hgtransformハンドル

function hEarthTransform = drawEarth3D(GMST0, const, sunDir)
arguments
    GMST0 (1,1) double
    const
    sunDir (1,3) double = []
end

libPath = '/Users/yyoshimula/Dropbox/MATLAB/yoshimuLibrary';

% 地球サーフェス
[xE, yE, zE] = ellipsoid(0, 0, 0, const.RE, const.RE, const.RE, 36);
earthSurf = surf(xE, yE, -zE, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
hEarthTransform = hgtransform;
set(earthSurf, 'Parent', hEarthTransform);
set(hEarthTransform, 'Matrix', makehgtform('zrotate', GMST0));

% テクスチャ
earthImg = imread(fullfile(libPath, 'visualization', 'naturalEarth.jpg'));
set(earthSurf, 'FaceColor', 'texturemap', 'CData', earthImg, 'EdgeColor', 'none');

% 影オーバーレイ（慣性座標系に固定）
if ~isempty(sunDir)
    [xS, yS, zS] = ellipsoid(0, 0, 0, const.RE*1.002, const.RE*1.002, const.RE*1.002, 72);
    shadowSurf = surf(xS, yS, -zS, 'EdgeColor', 'none');

    norms = sqrt(xS.^2 + yS.^2 + zS.^2);
    cosAngle = (xS./norms)*sunDir(1) + (yS./norms)*sunDir(2) + (-zS./norms)*sunDir(3);
    alphaData = max(0, -cosAngle) * 0.95;

    set(shadowSurf, 'FaceColor', 'k', 'FaceAlpha', 'flat', ...
        'AlphaData', alphaData, 'AlphaDataMapping', 'none');

    % 太陽方向の矢印
    arrowLen = const.RE * 2;
    quiver3(0, 0, 0, sunDir(1)*arrowLen, sunDir(2)*arrowLen, sunDir(3)*arrowLen, ...
        'Color', [1, 0.8, 0], 'LineWidth', 2, 'MaxHeadSize', 0.5);
    text(sunDir(1)*arrowLen*1.1, sunDir(2)*arrowLen*1.1, sunDir(3)*arrowLen*1.1, ...
        'Sun', 'Color', [1, 0.8, 0], 'FontSize', 10);
end

end
