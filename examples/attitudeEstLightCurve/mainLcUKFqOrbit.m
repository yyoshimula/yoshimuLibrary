%[text] # photometric attitude (quaternion) and orbit estimation using UKF, not yet
%[text] UKFでlight curve attitude and orbit estimation
%[text] ## state variables
%[text] ${\\bf x}=\[{\\bf r}^T, {\\bf q}^T\]^T$
%[text] `r`: position vetor at inertial frame
%[text] `q`: quaternion where q(4) is the scalar part 
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20220805 y.yoshimura, minor changes
%[text] 20210209  y.yoshimura, y.yoshimula@gmail.com
%[text] See also orbitConst, mainAttiOrbit.
clc
clear
cls
%[text] ## constants
earthVSOP = vsopConst; % VSOP parameters
const = orbitConst; % constants for orbital motion
%%
%[text] ## simulation conditions
observation = 'mag'; % light curve magnitude or its inverse value

% random number generator
rng(2022,'twister');
%[text] ## load true value
load('flatPlate.mat')
% ****.mat includes following variables
% t_:
% sat: satellite configuration
% mApp:
% w: angular rate, Nx3
% r: satellite position at inertial frame, m
% obsECI:
% jdHistory
tn_ = length(t_);
%     sat.nu = ones(length(sat.Cs),1) .* 300;
%     sat.nv = ones(length(sat.Cs),1) .* 300;
sigTrue = 0.1;
Rtrue = sigTrue^2; % observation noise variance for relative magnitude
%[text] ## UKF paras.
n_ = 6; % # of state variable := [position and error GRP]

% UKF Parameters
ukfPara.alp = 1e-4;
ukfPara.beta = 2;
ukfPara.kap = 3 - n_;
% ukaPara.kap = 0;
ukfPara.lam = ukfPara.alp^2 * (n_ + ukfPara.kap) - n_ %[output:58727e5b]
% ukfPara.lam = 0;

ukfPara.w0m = ukfPara.lam / (n_ + ukfPara.lam); % for mean
ukfPara.wim = 1 / (2 * (n_ + ukfPara.lam));
ukfPara.w0c = ukfPara.lam / (n_ + ukfPara.lam) + (1 - ukfPara.alp^2 + ukfPara.beta); % for covariance
ukfPara.wic = ukfPara.wim;

ukfPara.w0m + 2*n_*ukfPara.wim %[output:94fe9fa1]
ukfPara.w0c + 2*n_*ukfPara.wic %[output:557670d8]

%[text] ## oservation noise (standard deviation)
if strcmp(observation, 'mag')
    sigEst = 3 * sigTrue;
    Rest = sigEst.^2;
    Qest = 0 * deg2rad(1e-2)^2 .* eye(3); % process noise
else
    sigEst = 1e-8;
    Rest = sigEst.^2;
    Qest = 0 * deg2rad(1e-2)^2 .* eye(3); % process noise
end

% azimuth and elevation for observation
% true
azi = atan2(r(:,2), r(:,1));
ele = asin(r(:,3) ./ vecnorm(r,2,2));
% noise
aziObs = azi + deg2rad(1) .* randn(tn_,1);
eleObs = ele + deg2rad(1) .* randn(tn_,1);
%%
%[text] ### pre-allocation
xEst = zeros(tn_, n_); % [position, error GRP]
p = zeros(n_, n_, tn_); % covariance matrix (Pcov)
qGlobal = zeros(tn_, 4); % quaternion estimate
ye0 = 0;
yez = zeros(2*n_,1);
yeOut = zeros(tn_, 1);
 
% observations
mAppNoise = mApp + sigTrue .* randn(tn_,1); % true obs. (relative magnitude)
dTmp = vecnorm(obsRel, 2, 2); % m, distance between sat and observer
mAppInv = magInv(mAppNoise(:,1), dTmp); % true obs. (relative magnitude inverse)
%%
%[text] ## initial estimate
% initial estimate
qGlobal(1,:) = q(1,:);  % true initial
iniQerr = deg2rad(60);
% large initial error
qGlobal(1,:) = qMult(4, 1, zyx2q(4, iniQerr, iniQerr, iniQerr), q(1,:));
% small initial error
qGlobal(1,:) = qMult(4, 1, zyx2q(4, deg2rad(20), deg2rad(-10), deg2rad(-10)), q(1,:));

iniRerr = 10 * 10^3; % m
iniR = r(1,:); % true initial position, m
iniR = r(1,:) + iniRerr .* ones(1,3);
%[text] #### initial covariance for 3$\\sigma\n$bound 
p(:,:,1) = diag([(iniRerr/3*10^3)^2.*ones(1,3), (iniQerr/3)^2 .* ones(1,3)]);

% x:= [positiona, error GRP]
xEst(1,:) = [iniR, 0 0 0];
%[text] ### Gneralized Rodrigues parameters
aGRP = 1;
fGRP = 2 * (aGRP + 1);
%%
%[text] ## UKF
for i = 1:tn_-1 %[output:group:3c813d50]
    dt_ = t_(i+1) - t_(i);
    [xEst(i+1,:), qGlobal(i+1,:), p(:,:,i+1)] = lcAttiOrbitUKF(jdHistory(i+1), dt_, xEst(i,:), p(:,:,i), mAppNoise(i+1), ... %[output:8e7cb888]
        aziObs(i+1), eleObs(i+1), ukfPara, Qest, Rest, fGRP, aGRP, qGlobal(i,:), w(i,:), v(i,:), obsECI(i+1,:), sat, const, earthVSOP); %[output:8e7cb888]
end %[output:group:3c813d50]
%[text] ## data handling
%[text] ### orbit
rEst = xEst(:,1:3);
rErr = r - rEst;
%[text] ### attitude
qe = qErr(4, q, qGlobal);
zyxE = q2zyx(4, qe);
phiE = zyxE(:,1);
thetaE = zyxE(:,2);
psiE = zyxE(:,3);

%[text] ## show figures
%[text] ### orbit
figure
tiledlayout(3,1),nexttile
plot(t_, rErr(:,1), 'r')
nexttile
plot(t_, rErr(:,2), 'g')
nexttile
plot(t_, rErr(:,3), 'b')

figure
tiledlayout(3,1),nexttile
plot(t_, r(:,1), 'r-'), hold on
plot(t_, rEst(:,1), 'gx')
nexttile
plot(t_, r(:,2), 'r-'), hold on
plot(t_, rEst(:,2), 'gx')
nexttile
plot(t_, r(:,3), 'r-'), hold on
plot(t_, rEst(:,3), 'bx')
%[text] ### attitude
figure
tiledlayout(4,1),nexttile
plot(t_, qe(:,1), 'r')
nexttile
plot(t_, qe(:,2), 'g')
nexttile
plot(t_, qe(:,3), 'b')
nexttile
plot(t_, qe(:,4), 'k')

figure
tiledlayout(3,1),nexttile
plot(t_/60, rad2deg(phiE), 'r')
xlim([0 60*6]),ylim([-60 60])
nexttile
plot(t_/60, rad2deg(thetaE), 'g')
xlim([0 60*6]),ylim([-60 60])
nexttile
plot(t_/60, rad2deg(psiE), 'b')
xlim([0 60*6]),ylim([-60 60])

figure
plot(t_/60, mApp)
xlabel('time [min]'), ylabel('Light curve, relative magnitude')
set(gca, 'YDir', 'reverse')


%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":33.1}
%---
%[output:58727e5b]
%   data: {"dataType":"textualVariable","outputData":{"header":"フィールドをもつ struct:","name":"ukfPara","value":"     alp: 1.0000e-04\n    beta: 2\n     kap: -3\n     lam: -6.0000\n"}}
%---
%[output:94fe9fa1]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"1.0000"},"version":0}
%---
%[output:557670d8]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"4.0000"},"version":0}
%---
%[output:8e7cb888]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"次を使用中のエラー: <a href=\"matlab:matlab.internal.language.introspective.errorDocCallback('sun')\" style=\"font-weight:bold\">sun<\/a>\n出力引数が多すぎます。\n\nエラー: <a href=\"matlab:matlab.internal.language.introspective.errorDocCallback('lcAttiOrbitUKF', '\/Users\/yyoshimula\/Dropbox\/MATLAB\/MATLABdrive\/yoshimuLibrary\/examples\/estimationLightCurve\/lcAttiOrbitUKF.mlx', 50)\" style=\"font-weight:bold\">lcAttiOrbitUKF<\/a> (<a href=\"matlab: opentoline('\/Users\/yyoshimula\/Dropbox\/MATLAB\/MATLABdrive\/yoshimuLibrary\/examples\/estimationLightCurve\/lcAttiOrbitUKF.mlx',50,0)\">行 50<\/a>)\n[lon, lat, dist] = sun(jd, const, earthVSOP);"}}
%---
