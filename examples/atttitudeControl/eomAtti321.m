%[text] # Equations of motion using 3-2-1 Euler angles
function dxdt = eomAtti321(~, x_, para)
phi = x_(1);
theta = x_(2);
psi = x_(3);
wx = x_(4);
wy = x_(5);
wz = x_(6);

w = [wx; wy; wz];

MOI = para.MOI;
phiD = para.phiD;
wD = [0; 0; 0]; % desired angular rate

trq = [0; 0; 0]; % initialize
%[text] ## Control law
k1 = 5;
k2 = 5;

%[text] ### Euler angle feedback
tmp = phiD - [phi; theta; psi];
phiErr = [tmp(3); tmp(2); tmp(1)]; % yoshimuraの定義ではphiはz軸周りの角度なので入れ替え

wErr = wD - w;

trq = k1 .* phiErr + k2 .* wErr;

%[text] ### quaternion feedback
% q = zyx2q(4, phi, theta, psi);
% qd = zyx2q(4, phiD(1), phiD(2), phiD(3));
% qe = qErr(4, q, qd)';
% 
% trq = k1 .* qe(1:3) + k2 .* wErr;

%[text] ### dynamics
dwdt = -cross(w, MOI * w) + trq;
dwdt = MOI^(-1) * dwdt;

%[text] ### Kinematics
kine = [0, sin(psi), cos(psi)
    0, cos(psi)*cos(theta), -sin(psi)*cos(theta)
    cos(theta), sin(psi)*sin(theta), cos(psi)*sin(theta)];
kine = 1/cos(theta) .* kine * w;


dxdt = [kine
    dwdt];


end

%[appendix]{"version":"1.0"}
%---
