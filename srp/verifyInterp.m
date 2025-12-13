%% not yet completed
clc
clear
cls

format long
%% config
config.phiI = deg2rad(360 * rand());
config.phiI = deg2rad(-90);
config.dTheta = 0.1;
config.thetaIspan = deg2rad(0:config.dTheta:90);

load('correctionPara.mat')

%% satellite configuration
nFacet = 1;
sat.pos = [0 0 0];
sat.normal = [0 0 1];
sat.area = 1;
sat.qlb = [0 0 0 1];

sat.F0 = 0.5 .* ones(nFacet, 1);
sat.kappa = 0 .* ones(nFacet, 1);
sat.Cd = 0.5 .* ones(nFacet, 1);
sat.Cs = sat.F0;
sat.Ca = 0.0 .* ones(nFacet, 1);
sat.mCT = 0.2 .* ones(nFacet, 1);

%% constants & parameters
const = orbitConst;
d = au2km(1.0, const) * 10^3; % m, distance from sat to sun
dAU = km2au(d ./ 10^3, const); % AU
S0 = const.S0; % Solar constant, W/m^2
c = const.c; % light speed, m/s
coeff = -S0 / c / dAU^2;

%% verification
srpCsCT = zeros(length(config.thetaIspan), 3);
srpCsApprox = zeros(length(config.thetaIspan), 3);
for i = 1:length(config.thetaIspan)
    thetaI = config.thetaIspan(i);
    phiI = config.phiI;
    sunB = [sin(thetaI)*cos(phiI), sin(thetaI)*sin(phiI), cos(thetaI)]; % sun direction

    NS = sat.normal * sunB'; % nFacet x 1

    % true
    [sat, ~, srpCsCT(i,:)] = srpCT(sat, sunB, d, const);

    % corrected
    [~, ~, srpCsApprox(i,:)] = srpCTinterp(sat, sunB, d, const, correctionPara);

end

%% analytical approximate solutio for s = n
% lamCT = 2 / sat.mCT^2;
% srpCsAnalytic = 8 * pi / lamCT^4 * (lamCT^3 - 5 * lamCT^2 + 12 * lamCT - 12 + exp(lamCT/sqrt(2)*(1-sqrt(2))) * (2 * lamCT^2 - 6 * sqrt(2)*lamCT + 12));

% sunB = [0, 0, 1];
% rRef = 2 * (sat.normal * sunB') .* sat.normal - sunB;
% M = ctM(sat, rRef, sunB);
% srpCsAnalytic = srpCsAnalytic * sat.area * coeff * M;

% [srpCsCT(1,3), srpCsAnalytic]
% 100 * abs(srpCsCT(1,3) - srpCsAnalytic) / abs(srpCsCT(1,3))

%% show figs for comparison
figure
tiledlayout(3,1)
nexttile
plot(rad2deg(config.thetaIspan), srpCsCT(:,1), 'b', 'DisplayName', 'True'), hold on
plot(rad2deg(config.thetaIspan), srpCsApprox(:,1), 'r', 'DisplayName', 'Approximate')
legend
nexttile
plot(rad2deg(config.thetaIspan), srpCsCT(:,2), 'b', 'DisplayName', 'True'), hold on
plot(rad2deg(config.thetaIspan), srpCsApprox(:,2), 'r', 'DisplayName', 'Approximate')
legend
nexttile
plot(rad2deg(config.thetaIspan), srpCsCT(:,3), 'b', 'DisplayName', 'True'), hold on
plot(rad2deg(config.thetaIspan), srpCsApprox(:,3), 'r', 'DisplayName', 'Approximate')
legend

% error
absErr = abs(srpCsCT - srpCsApprox);
relErr = 100 .* absErr ./ abs(srpCsCT);
figure
tiledlayout(3,1)
nexttile
plot(rad2deg(config.thetaIspan), absErr(:,1), 'DisplayName', 'Absolute Error (x)')
nexttile
plot(rad2deg(config.thetaIspan), absErr(:,2), 'DisplayName', 'Absolute Error (y)')
nexttile
plot(rad2deg(config.thetaIspan), absErr(:,3), 'DisplayName', 'Absolute Error (z)')
sgtitle('absolute error')

figure
tiledlayout(3,1)
nexttile
plot(rad2deg(config.thetaIspan), relErr(:,1), 'DisplayName', 'Relative Error (x)')
ylim([0 100])
nexttile
plot(rad2deg(config.thetaIspan), relErr(:,2), 'DisplayName', 'Relative Error (y)')
ylim([0 100])
nexttile
plot(rad2deg(config.thetaIspan), relErr(:,3), 'DisplayName', 'Relative Error (z)')
ylim([0 100])
sgtitle('relative error')