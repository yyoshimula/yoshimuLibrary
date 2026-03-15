%% Example Orbits - 様々な軌道タイプの可視化例
% ISS, 静止軌道, Molniya軌道, 太陽同期軌道などの例
%
% Author: Generated with Claude Code
% Date: 2026-03-14

clear; close all; clc;

%% ライブラリパスの追加
libPath = '/Users/yyoshimula/Dropbox/MATLAB/yoshimuLibrary';
addpath(genpath(libPath));

const = orbitConst();

%% 軌道タイプの選択
fprintf('=== 軌道タイプを選択 ===\n');
fprintf('1: ISS (低軌道, 傾斜角51.6度)\n');
fprintf('2: 静止軌道 (GEO)\n');
fprintf('3: Molniya軌道 (高離心率楕円軌道)\n');
fprintf('4: 太陽同期軌道 (SSO)\n');
fprintf('5: GPS軌道 (中軌道, MEO)\n');
fprintf('6: カスタム軌道\n');
fprintf('7: 複数軌道を同時表示\n\n');

orbitType = input('選択 (1-7): ');

switch orbitType
    case 1
        %% ISS軌道
        oe = [const.RE + 420, 0.0001, deg2rad(51.6), deg2rad(30), 0, 0];
        visualizeOrbit(oe, 'nOrbits', 1.5, 'animSpeed', 100);

    case 2
        %% 静止軌道 (GEO)
        % 軌道長半径: 42164 km (高度約35786 km)
        a_geo = 42164;
        oe = [a_geo, 0.0001, deg2rad(0.1), deg2rad(0), 0, 0];
        visualizeOrbit(oe, 'nOrbits', 1, 'nPoints', 360, 'animSpeed', 30);

    case 3
        %% Molniya軌道
        % 周期12時間, 傾斜角63.4度, 高離心率
        T_molniya = 12 * 3600;  % 12時間
        a_molniya = (const.GE * (T_molniya / (2*pi))^2)^(1/3);
        oe = [a_molniya, 0.74, deg2rad(63.4), deg2rad(45), deg2rad(270), 0];
        visualizeOrbit(oe, 'nOrbits', 2, 'nPoints', 800, 'animSpeed', 80);

    case 4
        %% 太陽同期軌道 (SSO)
        % 高度約700km, 傾斜角約98度
        oe = [const.RE + 700, 0.001, deg2rad(98.2), deg2rad(0), deg2rad(90), 0];
        visualizeOrbit(oe, 'nOrbits', 2, 'animSpeed', 100);

    case 5
        %% GPS軌道 (MEO)
        % 周期約12時間, 傾斜角55度
        T_gps = 11.967 * 3600;
        a_gps = (const.GE * (T_gps / (2*pi))^2)^(1/3);
        oe = [a_gps, 0.01, deg2rad(55), deg2rad(0), 0, 0];
        visualizeOrbit(oe, 'nOrbits', 1, 'nPoints', 400, 'animSpeed', 50);

    case 6
        %% カスタム軌道
        fprintf('\n軌道要素を入力:\n');
        alt = input('高度 [km]: ');
        ecc = input('離心率 [-]: ');
        inc_deg = input('軌道傾斜角 [deg]: ');
        RAAN_deg = input('昇交点赤経 [deg]: ');

        oe = [const.RE + alt, ecc, deg2rad(inc_deg), deg2rad(RAAN_deg), 0, 0];
        visualizeOrbit(oe, 'nOrbits', 2);

    case 7
        %% 複数軌道の同時表示
        multiOrbitVisualization(const);

    otherwise
        error('無効な選択です');
end

%% 複数軌道同時表示関数
function multiOrbitVisualization(const)
% 複数の軌道を同一図にアニメーション表示

libPath = '/Users/yyoshimula/Dropbox/MATLAB/yoshimuLibrary';
addpath(genpath(libPath));

%% 軌道の定義
orbits = {
    struct('name', 'ISS',     'oe', [const.RE + 420, 0.0001, deg2rad(51.6), deg2rad(30),  0,            0], 'color', [0 1 1])
    struct('name', 'GEO',     'oe', [42164,          0.0001, deg2rad(0.1),   0,             0,            0], 'color', [1 1 0])
    struct('name', 'Molniya', 'oe', [26554,          0.74,   deg2rad(63.4),  deg2rad(120),  deg2rad(270), 0], 'color', [1 0 1])
    struct('name', 'GPS',     'oe', [26560,          0.01,   deg2rad(55),    deg2rad(240),  0,            0], 'color', [0 1 0])
    };

nOrbits = length(orbits);

%% 時刻設定（最長周期の軌道に合わせる）
nPoints = 600;
Tmax = 0;
for j = 1:nOrbits
    T = 2 * pi * sqrt(orbits{j}.oe(1)^3 / const.GE);
    Tmax = max(Tmax, T);
end
tSpan = linspace(0, Tmax, nPoints);

jd0 = gc2jd(2026, 8, 14, 12, 0, 0);
jdArray = jd0 + s2day(tSpan);
GMSTarray = gmst(jdArray');

%% 各軌道を伝播 + 太陽/影
allR = cell(nOrbits, 1);
allNu = cell(nOrbits, 1);
for j = 1:nOrbits
    [rI, ~] = propagateOrbit(orbits{j}.oe, tSpan, const.GE);
    allR{j} = rI;
    [~, allNu{j}] = calcSunShadow(jdArray, rI, const);
end

% 太陽方向（1回だけ計算）
[sunPos, ~] = calcSunShadow(jdArray, allR{1}, const);
sunDir = sunPos(1,:) / norm(sunPos(1,:));

%% 図の作成
fig = figure('Name', 'Multiple Orbits', 'Position', [50, 50, 1200, 900], 'Color', 'k');
hold on; axis equal; view(30, 20);

% 地球
hEarthTransform = drawEarth3D(GMSTarray(1), const, sunDir);

% 軌道線と衛星マーカー
hSats = gobjects(nOrbits, 1);
legends = {};
for j = 1:nOrbits
    rI = allR{j}; col = orbits{j}.color;
    sunlit = allNu{j} > 0.5;

    plot3(rI(sunlit,1), rI(sunlit,2), rI(sunlit,3), '.', 'Color', col, 'MarkerSize', 2);
    plot3(rI(~sunlit,1), rI(~sunlit,2), rI(~sunlit,3), '.', 'Color', col*0.3, 'MarkerSize', 2);

    hSats(j) = plot3(rI(1,1), rI(1,2), rI(1,3), 'o', ...
        'Color', col, 'MarkerSize', 10, 'MarkerFaceColor', col);
    legends{end+1} = orbits{j}.name;
end

xlabel('X [km]', 'Color', 'w'); ylabel('Y [km]', 'Color', 'w'); zlabel('Z [km]', 'Color', 'w');
title('様々な衛星軌道', 'Color', 'w', 'FontSize', 16);
ax = gca;
ax.Color = 'k'; ax.XColor = 'w'; ax.YColor = 'w'; ax.ZColor = 'w';
ax.GridColor = [0.3, 0.3, 0.3]; grid on;
legend(hSats, legends, 'TextColor', 'w', 'Location', 'northeast');

hTime = annotation('textbox', [0.35, 0.95, 0.3, 0.05], ...
    'String', '', 'Color', 'w', 'EdgeColor', 'none', ...
    'FontSize', 12, 'HorizontalAlignment', 'center');

%% アニメーション
fprintf('\nアニメーション再生中... (Ctrl+Cで停止)\n');
dt = [2026, 8, 14, 12, 0, 0];

for k = 1:nPoints
    if ~isvalid(fig), break; end

    for j = 1:nOrbits
        rI = allR{j};
        set(hSats(j), 'XData', rI(k,1), 'YData', rI(k,2), 'ZData', rI(k,3));
    end

    set(hEarthTransform, 'Matrix', makehgtform('zrotate', GMSTarray(k)));

    currentTime = datestr(datetime(dt(1),dt(2),dt(3),dt(4),dt(5),dt(6)) + ...
        seconds(tSpan(k)), 'yyyy/mm/dd HH:MM:SS');
    set(hTime, 'String', currentTime);

    drawnow;
    pause(0.02);
end

fprintf('アニメーション完了\n');
end
