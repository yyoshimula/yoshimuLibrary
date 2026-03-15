%% Satellite Orbit Visualization
% 衛星の軌道運動を可視化するプログラム
%
% 使い方:
%   visualizeOrbit()              % デフォルトのISS軌道
%   visualizeOrbit(oe)            % カスタム軌道要素を指定
%   visualizeOrbit(oe, options)   % オプション付き
%   visualizeOrbit([], rvFile='data.mat')  % 位置速度履歴ファイルから読み込み
%
% 軌道要素 oe = [a, e, i, RAAN, omega, f]
%   a: 軌道長半径 [km], e: 離心率, i: 傾斜角 [rad]
%   RAAN: 昇交点赤経 [rad], omega: 近点引数 [rad], f: 真近点角 [rad]
%
% rvFile (.mat) に必要な変数 (rI または oeArray のいずれか):
%   rI: 位置履歴 [km], Nx3 (ECI)
%   oeArray: 軌道要素履歴 [a,e,i,RAAN,omega,f], Nx6
%   tSpan: 時刻配列 [s], 1xN or Nx1
%   jd0: エポックのユリウス日 (スカラ, 省略時はstartDateを使用)
%   oeFlag: oeArrayの第6列の種別 (1=真近点角(default), 0=平均近点角)
%   vI: 速度履歴 [km/s], Nx3 (省略可)

function visualizeOrbit(oe, options)

arguments
    % デフォルト: ISS風の軌道
    oe double = [6798.137, 0.0001, 0.9006, 0.5236, 0, 0]
    options.nOrbits (1,1) double = 2
    options.nPoints (1,1) double = 500
    options.animate (1,1) logical = true
    options.showShadow (1,1) logical = true
    options.showGroundTrack (1,1) logical = true
    options.startDate (1,6) double = [2026 8 14 12 0 0]
    options.animSpeed (1,1) double = 50
    options.rvFile (1,:) char = ''
end

%% セットアップ
libPath = '/Users/yyoshimula/Dropbox/MATLAB/yoshimuLibrary';
addpath(genpath(libPath));
const = orbitConst();

%% 位置速度履歴ファイルからの読み込み or 軌道伝播
useFile = ~isempty(options.rvFile);

if useFile
    %% ファイルから読み込み
    if ~isfile(options.rvFile)
        error('ファイルが見つかりません: %s', options.rvFile);
    end
    data = load(options.rvFile);
    if ~isfield(data, 'tSpan')
        error('rvFileにはtSpan (1xN or Nx1)が必要です');
    end
    tSpan = data.tSpan(:)';
    if isfield(data, 'rI')
        rI = data.rI;
    elseif isfield(data, 'oeArray')
        oeFlag = 1;
        if isfield(data, 'oeFlag'), oeFlag = data.oeFlag; end
        [rI, ~] = oe2rv(data.oeArray, oeFlag, const.GE);
    else
        error('rvFileにはrI (Nx3)またはoeArray (Nx6)が必要です');
    end
    dt = options.startDate;
    if isfield(data, 'jd0')
        jd0 = data.jd0;
    else
        jd0 = gc2jd(dt(1), dt(2), dt(3), dt(4), dt(5), dt(6));
    end
    jdArray = jd0 + s2day(tSpan);

    % 軌道情報の概算表示
    alt = vecnorm(rI, 2, 2) - const.RE;
    fprintf('=== ファイル読み込み: %s ===\n', options.rvFile);
    fprintf('データ点数 %d | 高度 %.0f~%.0f km | 時間 %.1f min\n', ...
        size(rI,1), min(alt), max(alt), tSpan(end)/60);
else
    %% 軌道要素から伝播
    T = 2 * pi * sqrt(oe(1)^3 / const.GE);
    fprintf('=== 軌道パラメータ ===\n');
    fprintf('高度 %.2f km | 離心率 %.6f | 傾斜角 %.2f deg | 周期 %.2f min\n', ...
        oe(1) - const.RE, oe(2), rad2deg(oe(3)), T/60);

    dt = options.startDate;
    jd0 = gc2jd(dt(1), dt(2), dt(3), dt(4), dt(5), dt(6));
    tSpan = linspace(0, T * options.nOrbits, options.nPoints);
    jdArray = jd0 + s2day(tSpan);

    [rI, ~] = propagateOrbit(oe, tSpan, const.GE);
end

%% 太陽・影の計算
GMSTarray = gmst(jdArray');
sunDir = [];
if options.showShadow
    [sunPos, nu] = calcSunShadow(jdArray, rI, const);
    sunDir = sunPos(1,:) / norm(sunPos(1,:));
    sunlit = nu > 0.5;
end

%% 地上軌跡の計算
if options.showGroundTrack
    [lat, lon] = eci2latlon(rI, GMSTarray, const);
end

%% 図の作成
fig = figure('Name', 'Satellite Orbit Visualization', ...
    'Position', [100, 100, 1400, 600], 'Color', 'k');

%% サブプロット1: 3D軌道表示
ax1 = subplot(1, 2, 1);
hold on; axis equal; view(30, 20);
xlabel('X [km]', 'Color', 'w'); ylabel('Y [km]', 'Color', 'w'); zlabel('Z [km]', 'Color', 'w');
title('衛星軌道 (慣性座標系)', 'Color', 'w', 'FontSize', 14);
ax1.Color = 'k'; ax1.XColor = 'w'; ax1.YColor = 'w'; ax1.ZColor = 'w';
ax1.GridColor = [0.3, 0.3, 0.3]; grid on;

% 地球
hEarthTransform = drawEarth3D(GMSTarray(1), const, sunDir);

% 軌道線
if options.showShadow
    plot3(rI(sunlit,1), rI(sunlit,2), rI(sunlit,3), '.', 'Color', [1,0.8,0], 'MarkerSize', 3);
    plot3(rI(~sunlit,1), rI(~sunlit,2), rI(~sunlit,3), '.', 'Color', [0.3,0.3,0.5], 'MarkerSize', 3);
else
    plot3(rI(:,1), rI(:,2), rI(:,3), 'c-', 'LineWidth', 1.5);
end

% 衛星マーカー
hSat = plot3(rI(1,1), rI(1,2), rI(1,3), 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');

%% サブプロット2: 地上軌跡
ax2 = subplot(1, 2, 2);
if options.showGroundTrack
    setupGroundTrack();
    title('地上軌跡 (Ground Track)', 'Color', 'w', 'FontSize', 14);

    if options.showShadow
        scatter(rad2deg(lon(sunlit)), rad2deg(lat(sunlit)), 8, [1,0.8,0], 'filled');
        scatter(rad2deg(lon(~sunlit)), rad2deg(lat(~sunlit)), 8, [0.4,0.4,0.6], 'filled');
    else
        plot(rad2deg(lon), rad2deg(lat), 'c.', 'MarkerSize', 5);
    end

    hGT = plot(rad2deg(lon(1)), rad2deg(lat(1)), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
end

%% 情報テキスト
if useFile
    infoStr = sprintf('高度: %.0f~%.0f km | 時間: %.1f min', ...
        min(alt), max(alt), tSpan(end)/60);
else
    infoStr = sprintf('高度: %.0f km | 周期: %.1f min | 傾斜角: %.1f°', ...
        oe(1) - const.RE, T/60, rad2deg(oe(3)));
end
annotation('textbox', [0.02, 0.02, 0.3, 0.08], ...
    'String', infoStr, ...
    'Color', 'w', 'EdgeColor', 'none', 'FontSize', 10);

hTime = annotation('textbox', [0.4, 0.95, 0.2, 0.05], ...
    'String', '', 'Color', 'w', 'EdgeColor', 'none', ...
    'FontSize', 12, 'HorizontalAlignment', 'center');

%% アニメーション
if options.animate
    fprintf('\nアニメーション再生中... (Ctrl+Cで停止)\n');

    for k = 1:length(tSpan)
        if ~isvalid(fig), break; end

        set(hSat, 'XData', rI(k,1), 'YData', rI(k,2), 'ZData', rI(k,3));
        set(hEarthTransform, 'Matrix', makehgtform('zrotate', GMSTarray(k)));

        if options.showGroundTrack && isvalid(hGT)
            set(hGT, 'XData', rad2deg(lon(k)), 'YData', rad2deg(lat(k)));
        end

        currentTime = datestr(datetime(dt(1),dt(2),dt(3),dt(4),dt(5),dt(6)) + ...
            seconds(tSpan(k)), 'yyyy/mm/dd HH:MM:SS');
        set(hTime, 'String', currentTime);

        drawnow;
        pause(1/options.animSpeed);
    end
    fprintf('アニメーション完了\n');
end

end
