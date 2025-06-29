%[text] # Dual quaternion kinematics and dynamics under SRP
%[text] 並進と回転を別々にといたものと比較
clc
clear

%[text] ## parameters and initial condtioin
global m MOI mu iniJD

const = orbitConst;
mu = const.GE;

earthVSOP = vsopConst; % for earth orientation

% for time calculation
iau06 = readIAU06;
leapJD = leapS; % load database
fname = "EOP_20_C04_one_file_1962-now.txt";
eopDataAll = readEOP(fname);
%[text] ## satellite model
model = 1;

if model == 1 %[output:group:4b5c7254]
    objName = 'flatPlate.obj';    
    sat = readSC(objName); % satellite shape and optical parameters %[output:46ef68e0]

    sat.J = diag([10, 15, 20]); % moment of inertia, kgm^2
    sat.m = 50; % mass, kg
    sat.nu = ones(size(sat.faces,1),1) .* 800; %とりあえずベタうち
    sat.nv = ones(size(sat.faces,1),1) .* 800;
    sat.F0 = ones(size(sat.faces,1),1) .* 0.5;
    sat.kappa = 0.0;
elseif model == 2
    objName = 'boxWing.obj';
    nCompo = 3; % # of object in .obj file
    sat = readSC(objName, nCompo); % satellite shape and optical parameters

    sat.J = diag([45, 50, 60]); % moment of inertia, kgm^2
    sat.m = 100; % mass, kg
    sat.kappa = 0.0;
    %   optical parameters at wing (manually...)
    sat.Ca(193:end,:) = 0.1;
    sat.Cd(193:end,:) = 0.6;
    sat.Cs(193:end,:) = 0.3;
    sat.nu = ones(size(sat.faces,1),1) .* 800; %とりあえずベタうち
    sat.nv = ones(size(sat.faces,1),1) .* 800;
    sat.F0 = ones(size(sat.faces,1),1) .* 0.5;
elseif model == 3

    objName = 'oneBoxWing.obj';
    nCompo = 2; % # of object in .obj file
    sat = readSC(objName, nCompo); % satellite shape and optical parameters

    sat.J = diag([45, 50, 60]); % moment of inertia, kgm^2
    sat.m = 100; % mass, kg
    sat.kappa = 0.0;
    sat.nu = ones(size(sat.faces,1),1) .* 800; %とりあえずベタうち
    sat.nv = ones(size(sat.faces,1),1) .* 800;
    sat.F0 = ones(size(sat.faces,1),1) .* 0.5;

else
    error('invalid model flag');
end %[output:group:4b5c7254]
showSC(sat,"Normal","on") % visualize satellite
%[text] ### orbital motion
%[text] $\\oe = \[a,e,i,\\Omega,\\omega,f\]\n$
%[text] ### using TLE
tle = readTLE('JCSAT2.txt', const);
oe = tle.oe;
iniJD = tle.jd;

%[text] ### or manual setting
% oe = [const.RE+800, 0.0, deg2rad(0), 0, 0, 0];

[iniPos, iniVelI] = oe2rv(oe, 1, const.GE);
n = sqrt(const.GE / oe(1)^3);
T = 2 * pi / n;

iniQ = rand(1,4);
iniQ = iniQ ./ norm(iniQ);

iniW = [0.0, 0.001, 0.003];

iniDQ = pos2dq(1, 4, iniPos, iniQ) %[output:2a684843]
iniVelB = qRotation(4, iniVelI, iniQ); %@body-fixed frame
dualVel = [iniW, 0, iniVelB, 0];

option = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
tspan = 0:2:T*1;
t_ = tspan;
%[text] ## dual quaternion
tic
[t_, x_] = ode45(@(t,x) dqEOMsrp(t,x, const, sat, earthVSOP), tspan, [iniDQ dualVel], option);
x_ = ode4(@(t,x) dqEOMsrp(t,x, const, sat, earthVSOP), tspan, [iniDQ dualVel]');
toc
%[text] ### data handling
dq = x_(:,1:8);
dualVel = x_(:,9:16);

[r, q] = dq2pos(1, 4, dq);

figure, hold on;
plot(t_, dq(:,1));
plot(t_, dq(:,2));
plot(t_, dq(:,3));
plot(t_, dq(:,4));

figure
plot3(dualVel(:,1), dualVel(:,2), dualVel(:,3))

figure
plot3(r(:,1), r(:,2), r(:,3)), hold on
drawEarth(0, 0.7, const);
axis equal

%[text] ## decoupled motion
iniX = [iniQ, iniW, iniPos, iniVelI];
tic
% [t_, x_] = ode45(@(t,x) decEOMsrp(t,x, const, sat, earthVSOP), tspan, iniX, option);
x_ = ode4(@(t,x) decEOMsrp(t,x, const, sat, earthVSOP), tspan, iniX');
toc

decR = x_(:,8:10);

%[text] ### pseudo true motion
iniX = [iniQ, iniW, iniPos, iniVelI];
option = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
tic
[t_, true_] = ode89(@(t,x) decEOMsrp(t,x, const, sat, earthVSOP), tspan, iniX, option);

toc

trueR = true_(:,8:10);

%[text] ## compare
figure
plot3(decR(:,1), decR(:,2), decR(:,3)), hold on
drawEarth(0, 0.7, const);
axis equal

figure
tiledlayout(3,1),nexttile
plot(t_, trueR(:,1) - r(:,1), 'r'),nexttile
plot(t_, trueR(:,2) - r(:,2), 'g'),nexttile
plot(t_, trueR(:,3) - r(:,3), 'b')

figure
tiledlayout(3,1),nexttile
plot(t_, trueR(:,1) - decR(:,1), 'r'),nexttile
plot(t_, trueR(:,2) - decR(:,2), 'g'),nexttile
plot(t_, trueR(:,3) - decR(:,3), 'b')

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":49.8}
%---
%[output:46ef68e0]
%   data: {"dataType":"text","outputData":{"text":"number of vertices:\n     8\n\nnumber of faces:\n    12\n\n","truncated":false}}
%---
%[output:2a684843]
%   data: {"dataType":"matrix","outputData":{"columns":8,"exponent":"4","name":"iniDQ","rows":1,"type":"double","value":[["0.000053321544824","0.000059281724812","0.000008310956623","0.000059778072576","0.729110706539882","1.042924618075312","-0.431025284472530","-1.624700208638053"]]}}
%---
