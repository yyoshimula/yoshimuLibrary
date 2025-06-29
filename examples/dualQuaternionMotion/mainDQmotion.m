%[text] # Dual quaternion kinematics and dynamics
%[text] 並進と回転を別々にといたものと比較
clc
clear

%[text] ## parameters and initial condtioin
global m MOI mu
MOI = diag([100, 150, 200]);
m = 100;

const = orbitConst;
mu = const.GE;

% trivial motion
% iniPos = [0 0 0];
% iniVelI = [1 1 0]; % @inertial frame

%[text] ### orbital motion
%[text] $\\oe = \[a,e,i,\\Omega,\\omega,f\]\n$
%[text] ### using TLE
tle = readTLE('JCSAT2.txt', const);
oe = tle.oe;

%[text] ### or manual setting
% oe = [const.RE+800, 0.0, deg2rad(0), 0, 0, 0];

[iniPos, iniVelI] = oe2rv(oe, 1, const.GE);
n = sqrt(const.GE / oe(1)^3);
T = 2 * pi / n;

iniQ = rand(1,4);
iniQ = iniQ ./ norm(iniQ);

iniW = [0.0, 0.001, 0.003];

iniDQ = pos2dq(1, 4, iniPos, iniQ) %[output:788e4d6d]
iniVelB = qRotation(4, iniVelI, iniQ); %@body-fixed frame
dualVel = [iniW, 0, iniVelB, 0];

option = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
tspan = 0:2:T*1;
t_ = tspan;
%[text] ## dual quaternion
tic
[t_, x_] = ode45(@(t,x) dqEOM(t,x), tspan, [iniDQ dualVel], option);
% x_ = ode4(@dqEOM, tspan, [iniDQ dualVel]');
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
[t_, x_] = ode45(@(t,x) decEOM(t,x), tspan, iniX, option);
% x_ = ode4(@decEOM, tspan, iniX');
toc

decR = x_(:,8:10);

%[text] ### pseudo true motion
iniX = [iniQ, iniW, iniPos, iniVelI];
option = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
tic
[t_, true_] = ode89(@(t,x) decEOM(t,x), tspan, iniX, option);

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
%[output:788e4d6d]
%   data: {"dataType":"matrix","outputData":{"columns":8,"exponent":"4","name":"iniDQ","rows":1,"type":"double","value":[["0.000071325279661","0.000011001810611","0.000031412466038","0.000061684046737","1.181626881000892","0.857157012548428","-1.232405965678443","-0.891595788952103"]]}}
%---
