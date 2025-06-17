%[text] # calculating line of sight (LOS) from relative orbital elements (ROE) for angles-only navigation
%[text] `roe`: relative orbtial elements (ROE), nx6 matrix
%[text] $\\delta \\alpha =(a\_d-a)/a \\\\\n\\delta \\lambda = (u\_d - u ) + (\\Omega\_d - \\Omega) \\cos{i} \\\\\n\\delta e\_x = e\_{xd} - e\_{x} \\\\\n\\delta e\_y = e\_{yd} - e\_{y} \\\\\n\\delta i\_x = i\_d - i \\\\\n\\delta i\_y = (\\Omega\_d - \\Omega) \\sin{i}\n$
%[text] ## inputs
%[text] `chiefOE, deputyOE`: absolute orbital elements of chief and deputy, nx6 matrix
%[text] $a , e, i, \\Omega, w, f({\\rm or} \~M)$ (m or km, -, rad, rad, rad, rad)
%[text] `flag`: 1 = true anomaly, 0:= mean anomaly , scalar
%[text] `GE`: gravitational constant of the Earth, m or km (unit must be unified with the position and velocity)
%[text] ## outputs
%[text] `azi`: -T軸方向からR軸への角度, nx1 vector, rad
%[text] `ele`: R–T平面からN方向への角度, nx1 vector, rad 
%[text] ## note
%[text] 重力定数と位置・速度は単位を合わせること
%[text] ## references 
%[text] Sullivan, Joshua, Generalized Angles-Only Navigation Architecture for Autonomous Distributed Space Systems, Journal of Guidance, Control, and Dynamics.
%[text] ## revisions
%[text] 20211027  y.yoshimura, y.yoshimula@gmail.com
%[text] See also roe2losApprox.
function [azi, ele] = oe2los(chiefOE, deputyOE, flag, GE)
% arguments
%     chiefOE (:,6) {mustBeNumeric}
%     deputyOE (:,6) {mustBeNumeric}
%     flag (1,1) {mustBeMember(flag, [0, 1])}
%     GE (1,1) {mustBeNumeric}
% end
%[text] ### calculating relative distance from absolute orbital elements
[rC, ~] = oe2rv(chiefOE, flag, GE);
[rD, ~] = oe2rv(deputyOE, flag, GE);

if flag == 0
    fC = trueAnomaly(chiefOE(:,1), chiefOE(:,2), chiefOE(:,6));
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
qIjk2rtn = zxz2q(4, raanC, incC, uC);
rtn = qRotation(4, (rD - rC), qIjk2rtn); % realtive distance expressed with RTN frame
%[text] ## LOS
azi = atan2(rtn(:,1), -rtn(:,2)); % -T軸方向からR軸への角度
ele = atan2(rtn(:,3), vecnorm(rtn(:,1:2), 2 ,2)); % R–T平面からN方向への角度

end

%[appendix]{"version":"1.0"}
%---
