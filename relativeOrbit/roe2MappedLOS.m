%[text] # calculating approximated line of sight (LOS) from relative orbital elements (ROE) for angles-only navigation
%[text] LOS calculation using 1st order ROE mapping to RTN frame
%[text] `roe`: relative orbtial elements (ROE), nx6 matrix
%[text] $\\delta \\alpha =(a\_d-a)/a \\\\\n\\delta \\lambda = (u\_d - u ) + (\\Omega\_d - \\Omega) \\cos{i} \\\\\n\\delta e\_x = e\_{xd} - e\_{x} \\\\\n\\delta e\_y = e\_{yd} - e\_{y} \\\\\n\\delta i\_x = i\_d - i \\\\\n\\delta i\_y = (\\Omega\_d - \\Omega) \\sin{i}\n$
%[text] `chiefOE`: absolute orbital elements of chief, nx6 matrix
%[text] $a , e, i, \\Omega, w, f({\\rm or} \~M)$ (km or m, -, rad, rad, rad, rad)
%[text] `flag`: 1 = true anomaly, 0:= mean anomaly , scalar
%[text] `azi`: -T軸方向からR軸への角度, nx1 vector, rad
%[text] `ele`: R–T平面からN方向への角度, nx1 vector, rad 
%[text] `GE`: gravitational constant of the Earth, m or km (unit must be unified with the position and velocity)
%[text] ## note
%[text] 重力定数とsemi-major axisの単位を合わせること
%[text] ## references 
%[text] Di Mauro, G. 2019 Minimum-Fuel Control Strategy for Spacecraft Formation Reconfiguration via Finite-Time Maneuvers Journal of Guidance, Control, and Dynamics.
%[text] ## revisions
%[text] 20211027  y.yoshimura, y.yoshimula@gmail.com
%[text] See also oe2roe.
function [azi, ele] = roe2mappedLOS(roe, chiefOE, flag, GE, chiefQbo)
% arguments
%     roe (:,6) {mustBeNumeric}
%     chiefOE (:,6) {mustBeNumeric}
%     flag (1,1) {mustBeMember(flag, [1, 0])}
%     GE (1,1) {mustBeNumeric}
%     chiefQbo (1,4) {mustBeNumeric}
% end

%[text] ### mapping from ROE to RTN
rel = roe2rtn(roe, chiefOE, flag, GE); % 1x6 vector

if (nargin < 5) % if chiefQbo is not provided, LOS is calculated in RTN frame
% do nothing

else % if chiefQbo is provided, LOS is calculated in body-fixed frame of chief

    rel = qRotation(4, rel, chiefQbo); % realtive distance expressed with body-fixed frame of chief

end

azi = atan2(rel(:,2), rel(:,1)); % R-T平面内の角度
ele = asin(rel(:,3) ./ vecnorm(rel, 2, 2)); % R-T平面からN方向への角度

end

%[appendix]{"version":"1.0"}
%---
