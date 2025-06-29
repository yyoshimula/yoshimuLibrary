%[text] # dual quaternion kinematics and dynamics
%[text] `scalar:` specify the definition of the quaternion 
%[text] `scalar == 0`
%[text] ${\\bf q} = \[q\_0,q\_1,q\_2,q\_3\]^T=\[\\cos(\\theta/2), {\\bf e}^T\\sin(\\theta/2)\]^T$
%[text] `scalar == 4`
%[text] ${\\bf q} = \[q\_1,q\_2,q\_3, q\_4\]^T=\[{\\bf e}^T\\sin(\\theta/2), \\cos(\\theta/2)\]^T$
%[text] x: dual quaternion \[real, dual\], 1x8
%[text] ## note
%[text] NA
%[text] ## references 
%[text] Sveier, A., & Egeland, O. (2020). Dual Quaternion Particle Filtering for Pose Estimation. IEEE Transactions on Control Systems Technology, 1-14.
%[text] ## revisions
%[text] 20231219  y.yoshimura y.yoshimula@gmail.com
%[text] See also pos2dq.
function dxdt = dqEOMsrp(t_, x_, const, sat, earthVSOP)
global iniJD

jd = iniJD + s2day(t_);

mu = const.GE;
m = sat.m;
MOI = sat.J;

%[text] ### kinematics
dq = x_(1:8); % column vector
w = x_(9:11);
dualVel = x_(9:16);
 
dqKine = 0.5 .* dqMult(4, 1, dualVel', dq')';
%[text] ### dynamics
J12 = [m.*eye(3), zeros(3,1)
    zeros(1,3), 1];

J21 = [MOI, zeros(3,1)
    zeros(1,3), 1];

J = [zeros(4,4), J12
    J21, zeros(4,4)];


% gravity
[r, q] = dq2pos(1, 4, dq');

aI = -mu / norm(r(1:3))^3 .* r(1:3);
fI = aI .* m; % 加速度でなくforceで与える

%[text] ## SRP
%[text] ### sun
sunPosI = sun(jd, const, earthVSOP); % km
sunPosB = qRotation(4, sunPosI, q);

sunD = norm(sunPosI) * 10^3; %m

% nu = shadow(r, sunPosI, const.RS, const.RE); % calc Earth's shadow

sunB = sunPosB ./ norm(sunPosB); % 1x3 sun directional (unit) vector

% simple model
[sat, ~, ~] = srpSimple(sat, sunB, sunD, const);
frc = sum(sat.force, 1)';
trq = sum(sat.torque, 1)';


%[text] ### dynamics (cont.)

fb = qRotation(4, fI, q);
fb = fb(:);

tb = zeros(3,1);

% add SRP force and torque
fb = fb + frc;
tb = tb + trq;

dualForce = [fb; 0; tb; 0];

tmp = dqMult(4, 0, dualVel', (J * dualVel)');
dqDyna = J^(-1) * (-tmp' + dualForce);
dqDyna(4) = 0;
dqDyna(8) = 0;

dxdt = [dqKine
    dqDyna];


end

%[appendix]{"version":"1.0"}
%---
