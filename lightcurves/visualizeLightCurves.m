%[text] # BRDF visualization
%[text] 20220406  y.yoshimura
%[text] (c) 2022 yasuhiro yoshimura
clc
clear
cls
%[text] ## object parameters
sat.pos = [0 0 0]; % position vector
sat.normal = [0 0 1]; % normal vector
sat.area = 1; % area, m^2

sat.F0 = 0.5;
sat.rho = 0.5;
sat.Cd = 0.5;
%[text] ### for Ashikhmin–Shirley model
sat.nu = 800; % 100 or 800
sat.nv = 800;
sat.uu = [1 0 0];
sat.uv = [0 1 0];
%[text] ### for Cook–Torrance model
sat.mCT = 0.1; % roughness parameter
%[text] ## illumination condition
const = orbitConst;

% sun direction
phiI = deg2rad(-90);
thetaI = deg2rad(0);
s = [cos(phiI) * sin(thetaI), sin(phiI) * sin(thetaI), cos(thetaI)];
sunPosI = au2km(1, const) * 10^3 .* s; % m 1x3 vector

% reference vector
azi = 0:deg2rad(2):2*pi;
ele = 0:deg2rad(2):deg2rad(80); % 90 degまですると値が大きくなりすぎるので．
[phiR, thetaR] = meshgrid(azi, ele);

x = cos(phiR) .* sin(thetaR);
y = sin(phiR) .* sin(thetaR);
z = cos(thetaR);

% perfect specular direction
rRef = 2 * (s * sat.normal') * sat.normal - s; % 1x3

% object attitude
q = [0 0 0 1];
satPos = [0 0 0];
%[text] ## true AS model
mApp = zeros(size(phiR));
for i = 1:size(phiR,1) %[output:group:89454a1f]
    for j = 1:size(phiR, 2)
        phi = phiR(i,j);
        theta = thetaR(i,j);
        % observer direction
        v = [cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)];
        obsPos = const.REm .* v;

        [mApp(i,j), ~] = lc(sat, 4, q, satPos, obsPos, sunPosI, 1, 'AS');         %[output:747dd0ea]

    end
end %[output:group:89454a1f]
%%
%[text] ## show figs
%[text] on hemisphere
l = 1.5; % arrow scale

figure
figSurf = surf(x, y, z, mApp, 'edgecolor', 'none')
colorbar
% alpha 0.4
axis equal
hold on
quiver3(0,0,0, l*s(1), l*s(2), l*s(3), 'r')
quiver3(0,0,0, l*sat.normal(1), l*sat.normal(2), l*sat.normal(3), 'g')
xlabel('x'), ylabel('y'), zlabel('z')
view([270+rad2deg(phiI) 30])
view([124.61 34.81])

% figure
% contourf(cos(phiR).*sin(thetaR), sin(phiR).*sin(thetaR), D);
% hold on
% plot(vNominal(1), vNominal(2), 'yx')
% colorbar
% axis equal
% xlabel('x'), ylabel('y')

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":44}
%---
%[output:747dd0ea]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"次を使用中のエラー: <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback('lc', '\/Users\/yyoshimula\/Dropbox\/MATLAB\/yoshimuLibrary\/lightcurves\/lc.m', 28)\" style=\"font-weight:bold\">lc<\/a> (<a href=\"matlab: opentoline('\/Users\/yyoshimula\/Dropbox\/MATLAB\/yoshimuLibrary\/lightcurves\/lc.m',28,0)\">行 28<\/a>)\n位置 8 の引数が無効です。 名前と値の引数に名前と値の両方が含まれていることを確認してください。"}}
%---
