% calculating line of sight (LOS) from relative orbital elements (ROE) for angles-only navigation
% roe: relative orbtial elements (ROE), nx6 matrix
% 
% inputs
% chiefOE, deputyOE: absolute orbital elements of chief and deputy, nx6 matrix
%  (m or km, -, rad, rad, rad, rad)
% flag: 1 = true anomaly, 0:= mean anomaly , scalar
% GE: gravitational constant of the Earth, m or km (unit must be unified with the position and velocity)
% outputs
% azi: -T軸方向からR軸への角度, nx1 vector, rad
% ele: R–T平面からN方向への角度, nx1 vector, rad 
% note
% 重力定数と位置・速度は単位を合わせること
% references 
% Sullivan, Joshua, Generalized Angles-Only Navigation Architecture for Autonomous Distributed Space Systems, Journal of Guidance, Control, and Dynamics.
% revisions
% 20211027  y.yoshimura, y.yoshimula@gmail.com
% See also roe2losApprox.

function [azi, ele] = oe2los_(chiefOE, deputyOE, flag, GE)

% calculating relative distance from absolute orbital elements
[rC, ~] = oe2rv_(chiefOE, flag, GE);
[rD, ~] = oe2rv_(deputyOE, flag, GE);

if flag == 0
    fC = trueAnomaly_(chiefOE(:,1), chiefOE(:,2), chiefOE(:,6));
else
    fC = chiefOE(:,6);
end
fC = mod(fC, 2*pi);

incC = chiefOE(:,3);
raanC = chiefOE(:,4);
raanC = mod(raanC, 2*pi);
uC = chiefOE(:,5) + fC;
uC = mod(uC, 2*pi);

% DCM from inertial frame to RTN frame
qIjk2rtn = zxz2q_(4, raanC, incC, uC);
rtn = qRotation_(4, (rD - rC), qIjk2rtn); % realtive distance expressed with RTN frame

% LOS
azi = atan2(rtn(:,1), -rtn(:,2)); % -T軸方向からR軸への角度
ele = atan2(rtn(:,3), vecnorm(rtn(:,1:2), 2 ,2)); % R–T平面からN方向への角度

end
