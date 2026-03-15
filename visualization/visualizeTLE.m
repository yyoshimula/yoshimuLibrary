%% Visualize Orbit from TLE or RV history files
% TLEファイルまたは位置速度履歴ファイルから衛星軌道を読み込み可視化
%
% 使い方:
%   visualizeTLE()                              % デフォルト(sample_tle.txt)
%   visualizeTLE('iss.tle')                     % 1ファイル
%   visualizeTLE('iss.tle', 'hubble.tle')       % 複数ファイル
%   visualizeTLE('iss.tle', 'nOrbits', 3)       % オプション付き
%   visualizeTLE('orbit.mat')                   % 位置速度履歴ファイル
%   visualizeTLE('iss.tle', 'orbit.mat')        % TLEと履歴ファイルの混在
%
% .matファイルに必要な変数 (rI または oeArray のいずれか):
%   rI: 位置履歴 [km], Nx3 (ECI)
%   oeArray: 軌道要素履歴 [a,e,i,RAAN,omega,f], Nx6
%   tSpan: 時刻配列 [s], 1xN or Nx1
%   jd0: エポックのユリウス日 (スカラ, 省略時は最初のTLEのエポックを使用)
%   oeFlag: oeArrayの第6列の種別 (1=真近点角(default), 0=平均近点角)
%   name: 衛星名 (文字列, 省略時はファイル名を使用)

function visualizeTLE(varargin)

%% 引数の解析（TLEファイル名とオプションを分離）
tleFiles = {};
optArgs = {};
optNames = {'nOrbits','nPoints','animate','showShadow','animSpeed'};
k = 1;
while k <= length(varargin)
    arg = char(varargin{k});
    if any(strcmp(arg, optNames))
        optArgs = varargin(k:end);
        break;
    elseif endsWith(arg, {'.txt','.tle','.TLE','.mat'})
        tleFiles{end+1} = arg;
    end
    k = k + 1;
end
if isempty(tleFiles), tleFiles = {'sample_tle.txt'}; end

p = inputParser;
addParameter(p, 'nOrbits', 2);
addParameter(p, 'nPoints', 500);
addParameter(p, 'animate', true);
addParameter(p, 'showShadow', true);
addParameter(p, 'animSpeed', 50);
parse(p, optArgs{:});
opts = p.Results;

%% セットアップ
libPath = '/Users/yyoshimula/Dropbox/MATLAB/yoshimuLibrary';
addpath(genpath(libPath));
const = orbitConst();
nFiles = length(tleFiles);

% 色の定義
colors = [0 1 1; 1 0.8 0; 0 1 0; 1 0 1; 1 0.5 0; 0.5 0.5 1];
if nFiles > size(colors,1)
    colors = [colors; rand(nFiles - size(colors,1), 3)];
end

%% TLE/.mat読み込み
sats = struct('name', {}, 'oe', {}, 'jd0', {}, 'T', {}, 'color', {}, ...
    'isFile', {}, 'rI', {});
for j = 1:nFiles
    if ~isfile(tleFiles{j}), error('ファイルが見つかりません: %s', tleFiles{j}); end

    if endsWith(tleFiles{j}, '.mat')
        % 位置速度履歴ファイルからの読み込み
        data = load(tleFiles{j});
        if ~isfield(data, 'tSpan')
            error('%s にはtSpan (1xN or Nx1)が必要です', tleFiles{j});
        end
        if ~isfield(data, 'rI') && ~isfield(data, 'oeArray')
            error('%s にはrI (Nx3)またはoeArray (Nx6)が必要です', tleFiles{j});
        end
        % oeArrayからrIを計算
        if ~isfield(data, 'rI')
            oeFlag = 1;
            if isfield(data, 'oeFlag'), oeFlag = data.oeFlag; end
            [data.rI, ~] = oe2rv(data.oeArray, oeFlag, const.GE);
        end
        if isfield(data, 'name')
            name = data.name;
        else
            [~, name, ~] = fileparts(tleFiles{j});
        end
        tSpanFile = data.tSpan(:)';
        T = tSpanFile(end);

        alt = vecnorm(data.rI, 2, 2) - const.RE;
        fprintf('%s (file): 高度 %.0f~%.0f km, 時間 %.1f min\n', ...
            name, min(alt), max(alt), T/60);

        sats(j).name = name;
        sats(j).oe = [];
        sats(j).jd0 = [];
        if isfield(data, 'jd0'), sats(j).jd0 = data.jd0; end
        sats(j).T = T;
        sats(j).color = colors(j,:);
        sats(j).isFile = true;
        sats(j).rI = data.rI;
        sats(j).tSpanFile = tSpanFile;
    else
        % TLEからの読み込み
        tle = readTLE(tleFiles{j}, const, 'yoshimuLibrary');
        oe = tle.oe(1,:);

        if isfield(tle, 'satName')
            name = strtrim(tle.satName(1,:));
        else
            [~, name, ~] = fileparts(tleFiles{j});
        end

        T = 2*pi*sqrt(oe(1)^3 / const.GE);
        fprintf('%s: 高度 %.0f km, 周期 %.1f min, 傾斜角 %.1f deg\n', ...
            name, oe(1)-const.RE, T/60, rad2deg(oe(3)));

        sats(j).name = name; sats(j).oe = oe;
        sats(j).jd0 = tle.jd(1); sats(j).T = T; sats(j).color = colors(j,:);
        sats(j).isFile = false;
        sats(j).rI = [];
        sats(j).tSpanFile = [];
    end
end

%% 共通の時間軸
Tmax = max([sats.T]);
nPoints = opts.nPoints;
tSpan = linspace(0, Tmax * opts.nOrbits, nPoints);

% jd0の決定: 最初に見つかったjd0を使用
jd0 = [];
for j = 1:nFiles
    if ~isempty(sats(j).jd0)
        jd0 = sats(j).jd0;
        break;
    end
end
if isempty(jd0)
    % jd0が見つからない場合はstartDateオプションのデフォルトを使用
    jd0 = gc2jd(2026, 8, 14, 12, 0, 0);
end
jdArray = jd0 + s2day(tSpan);
GMSTarray = gmst(jdArray');

%% 太陽位置
[sunPos, ~] = calcSunShadow(jdArray, zeros(nPoints,3), const);
sunDir = sunPos(1,:) / norm(sunPos(1,:));

%% 各衛星の軌道伝播 / ファイルデータの補間
for j = 1:nFiles
    if sats(j).isFile
        % ファイルデータを共通時間軸に補間
        rIorig = sats(j).rI;
        tOrig = sats(j).tSpanFile;
        rI = interp1(tOrig, rIorig, tSpan, 'spline', 'extrap');
        sats(j).rI = rI;
    else
        [rI, ~] = propagateOrbit(sats(j).oe, tSpan, const.GE);
        sats(j).rI = rI;
    end
    [~, sats(j).nu] = calcSunShadow(jdArray, sats(j).rI, const);
    [sats(j).lat, sats(j).lon] = eci2latlon(sats(j).rI, GMSTarray, const);
end

%% 図の作成
fig = figure('Name', 'TLE Orbit Visualization', ...
    'Position', [100, 100, 1400, 600], 'Color', 'k');

%% サブプロット1: 3D軌道表示
ax1 = subplot(1, 2, 1);
hold on; axis equal; view(30, 20);
xlabel('X [km]', 'Color', 'w'); ylabel('Y [km]', 'Color', 'w'); zlabel('Z [km]', 'Color', 'w');
title('衛星軌道 (慣性座標系)', 'Color', 'w', 'FontSize', 14);
ax1.Color = 'k'; ax1.XColor = 'w'; ax1.YColor = 'w'; ax1.ZColor = 'w';
ax1.GridColor = [0.3, 0.3, 0.3]; grid on;

hEarthTransform = drawEarth3D(GMSTarray(1), const, sunDir);

% 軌道線と衛星マーカー
hSats3D = gobjects(nFiles, 1);
legends = {};
for j = 1:nFiles
    rI = sats(j).rI; col = sats(j).color;
    sunlit = sats(j).nu > 0.5;

    plot3(rI(sunlit,1), rI(sunlit,2), rI(sunlit,3), '.', 'Color', col, 'MarkerSize', 2);
    plot3(rI(~sunlit,1), rI(~sunlit,2), rI(~sunlit,3), '.', 'Color', col*0.3, 'MarkerSize', 2);

    hSats3D(j) = plot3(rI(1,1), rI(1,2), rI(1,3), 'o', ...
        'Color', col, 'MarkerSize', 10, 'MarkerFaceColor', col);
    legends{end+1} = sats(j).name;
end
legend(hSats3D, legends, 'TextColor', 'w', 'Location', 'northeast');

%% サブプロット2: 地上軌跡
ax2 = subplot(1, 2, 2);
setupGroundTrack();
title('地上軌跡 (Ground Track)', 'Color', 'w', 'FontSize', 14);

hSatsGT = gobjects(nFiles, 1);
for j = 1:nFiles
    col = sats(j).color;
    sunlit = sats(j).nu > 0.5;

    scatter(rad2deg(sats(j).lon(sunlit)), rad2deg(sats(j).lat(sunlit)), 5, col, 'filled');
    scatter(rad2deg(sats(j).lon(~sunlit)), rad2deg(sats(j).lat(~sunlit)), 5, col*0.3, 'filled');

    hSatsGT(j) = plot(rad2deg(sats(j).lon(1)), rad2deg(sats(j).lat(1)), 'o', ...
        'Color', col, 'MarkerSize', 8, 'MarkerFaceColor', col);
end

%% 時刻表示
hTime = annotation('textbox', [0.4, 0.95, 0.2, 0.05], ...
    'String', '', 'Color', 'w', 'EdgeColor', 'none', ...
    'FontSize', 12, 'HorizontalAlignment', 'center');

%% アニメーション
if opts.animate
    fprintf('\nアニメーション再生中... (Ctrl+Cで停止)\n');

    for k = 1:nPoints
        if ~isvalid(fig), break; end

        for j = 1:nFiles
            rI = sats(j).rI;
            set(hSats3D(j), 'XData', rI(k,1), 'YData', rI(k,2), 'ZData', rI(k,3));
            if isvalid(hSatsGT(j))
                set(hSatsGT(j), 'XData', rad2deg(sats(j).lon(k)), ...
                    'YData', rad2deg(sats(j).lat(k)));
            end
        end

        set(hEarthTransform, 'Matrix', makehgtform('zrotate', GMSTarray(k)));

        [yr, mo, dy, hr, mi, sc] = jd2gc(jdArray(k));
        set(hTime, 'String', sprintf('%04d/%02d/%02d %02d:%02d:%02.0f', yr, mo, dy, hr, mi, sc));

        drawnow;
        pause(1/opts.animSpeed);
    end
    fprintf('アニメーション完了\n');
end

end
