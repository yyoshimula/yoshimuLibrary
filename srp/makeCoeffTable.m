clc
clear
cls

format long

%% configuration
config.gFigs = 1;
config.mSpan = [0.05, 0.145, 0.275, 0.45];
config.dTheta = 3;
config.thetaIspan = deg2rad(0:config.dTheta:90);
[config.thetaICases, config.mCases] = meshgrid(config.thetaIspan, config.mSpan); % mSpan x thetaIspan

% sun direction
config.phiI = deg2rad(-90);

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

%% constants & parameters
const = orbitConst;
d = au2km(1.0, const) * 10^3; % m, distance from sat to sun
dAU = km2au(d ./ 10^3, const); % AU
S0 = const.S0; % Solar constant, W/m^2
c = const.c; % light speed, m/s
coeff = -S0 / c / dAU^2;

%% specular

nMC = 10^6; % for Monte-Carlo integration

for i = 1:size(config.mCases, 1) % for each m
    for j = 1:size(config.mCases, 2) % for each thetaI
        sat.mCT = config.mCases(i,j); % roughness

        sunB = [sin(config.thetaICases(i,j))*cos(config.phiI), sin(config.thetaICases(i,j))*sin(config.phiI), cos(config.thetaICases(i,j))]; % sun direction

        NS = sat.normal * sunB';
        coeff = -S0 / c / dAU^2 * sat.area * NS;

        %% true (numerically)
        [tmp, srpCdCT(j,:,i), srpCsCT(j,:,i)] = srpCT(sat, sunB, d, const, 'Beckmann', nMC);
        ctSRP(j,:,i) = tmp.force;

        if i == 1 %% Perfect spceular (mは関係ないので)
            [tmp, srpCdSimple(j,:), srpCsSimple(j,:)] = srpSimple(sat, sunB, d, const);
            simpleSRP(j,:) = tmp.force;
        end

        % calc. correction parameters
        deltaS(i,j) = (srpCsCT(j,:,i) * sunB' - srpCsCT(j,:,i) * sat.normal' * NS) / (sat.Cs * NS^2 - sat.Cs);
        deltaS(i,j) = deltaS(i,j) / coeff;

        deltaN(i,j) = (srpCsCT(j,:,i) * sunB' * NS - srpCsCT(j,:,i) * sat.normal') / (2 * sat.Cs * NS^3 - 2 * sat.Cs * NS);
        deltaN(i,j) = deltaN(i,j) / coeff;

    end
end

%% gradient of correction parameters
nTheta = size(config.thetaICases, 2);
for i = 1:size(config.mCases, 1)
    gradDeltaS(i,2:nTheta-1) = (deltaS(i,3:nTheta) - deltaS(i,1:nTheta-2)) / (2 * deg2rad(config.dTheta)); % 中心差分
    gradDeltaN(i,2:nTheta-1) = (deltaN(i,3:nTheta) - deltaN(i,1:nTheta-2)) / (2 * deg2rad(config.dTheta)); % 中心差分

    gradDeltaS(i,1) = (-3*deltaS(i,1) + 4 * deltaS(i,2) - deltaS(i,3)) / (1 * deg2rad(config.dTheta)); % 2次精度の片側差分
    gradDeltaN(i,1) = (-3*deltaN(i,1) + 4 * deltaN(i,2) - deltaN(i,3)) / (1 * deg2rad(config.dTheta)); % 2次精度の片側差分

    gradDeltaS(i,nTheta) = (3 * deltaS(i,nTheta) - 4 * deltaS(i,nTheta-1) + deltaS(i,nTheta-2)) / deg2rad(config.dTheta); % 2次精度の片側差分
    gradDeltaN(i,nTheta) = (3 * deltaN(i,nTheta) - 4 * deltaN(i,nTheta-1) + deltaN(i,nTheta-2)) / deg2rad(config.dTheta); % 2次精度の片側差分
end


%% show figures
% comparison a single roughness
% for i = 1:length(mSpan)
%     figure
%     subplot(3,1,1)
%     plot(rad2deg(thetaICases), ctSRP(:,1,i), 'r'); hold on;
%     plot(rad2deg(thetaICases), simpleSRP(:,1), 'r--');
%     ylabel('x');
%     subplot(3,1,2)
%     plot(rad2deg(thetaICases), ctSRP(:,2,i), 'g');
%     plot(rad2deg(thetaICases), simpleSRP(:,2), 'g--');
%     ylabel('y');
%     subplot(3,1,3)
%     plot(rad2deg(thetaICases), ctSRP(:,3,i), 'b'); hold on;
%     plot(rad2deg(thetaICases), simpleSRP(:,3), 'b--');
%     ylabel('z');
%     xlabel('thetaI [deg]');
%     legend('Cook-Torrance', 'Simple', 'Location', 'best');
%     sgtitle(sprintf('Cook-Torrance SRP (m = %f)', mSpan(i)));
% end

% SRP norm comparison
% figure
% hold on
% plot(rad2deg(thetaICases), vecnorm(simpleSRP, 2, 2), '--', ...
%     'DisplayName', 'Simple');
% for i = 1:length(mSpan)
%     plot(rad2deg(thetaICases), vecnorm(ctSRP(:,:,i), 2, 2), ...
%         'DisplayName', sprintf('CT (m = %.3f)', mSpan(i)));
% end
% legend show

%% verify correction parameters
fCorrected = zeros(size(config.thetaICases, 2), 3, size(config.mCases, 1));
for i = 1:size(config.mCases, 1) % for each m
    for j = 1:size(config.thetaICases, 2) % for each thetaI
        sunB = [sin(config.thetaICases(i,j))*cos(config.phiI), sin(config.thetaICases(i,j))*sin(config.phiI), cos(config.thetaICases(i,j))]; % sun direction
        NS = sat.normal * sunB';
        coeff = -S0 / c / dAU^2 * sat.area * NS;
        fCorrected(j,:,i) = (1 - deltaS(i,j) * sat.Cs) .* sunB + (2/3*sat.Cd + 2 * NS * deltaN(i,j) * sat.Cs) .* sat.normal;
        fCorrected(j,:,i) = coeff * fCorrected(j,:,i);
    end
end

for i = 1:size(config.mCases, 1)
    figure
    tiledlayout(3,1), nexttile
    plot(rad2deg(config.thetaICases), ctSRP(:,1,i), 'r'); hold on
    plot(rad2deg(config.thetaICases), fCorrected(:,1,i), 'rx');
    nexttile
    plot(rad2deg(config.thetaICases), ctSRP(:,2,i), 'g'); hold on
    plot(rad2deg(config.thetaICases), fCorrected(:,2,i), 'gx');
    nexttile
    plot(rad2deg(config.thetaICases), ctSRP(:,3,i), 'b'); hold on
    plot(rad2deg(config.thetaICases), fCorrected(:,3,i), 'bx');
    legend('True', 'Corrected', 'Location', 'best');
    sgtitle(sprintf('Corrected SRP (m = %f)', config.mCases(i)));

    %% relative error
    % relErr(:,:,i) = 100 .* abs(ctSRP(:,:,i) - fCorrected(:,:,i)) ./ abs(ctSRP(:,:,i));
    % figure
    % tiledlayout(3,1), nexttile
    % plot(rad2deg(config.thetaICases), relErr(:,1,i), 'r');
    % nexttile
    % plot(rad2deg(config.thetaICases), relErr(:,2,i), 'g');
    % nexttile
    % plot(rad2deg(config.thetaICases), relErr(:,3,i), 'b');
    % sgtitle(sprintf('Relative Error (m = %f)', config.mCases(i)));
end

correctionPara = struct('deltaS', deltaS, 'deltaN', deltaN, ...
    'mSpan', config.mSpan, 'thetaIspan', config.thetaIspan);
save('correctionPara.mat', 'correctionPara');

% gFigs(config.gFigs)