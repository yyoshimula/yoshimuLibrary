%[text] # calculating approximated line of sight (LOS) from relative orbital elements (ROE) for angles-only navigation
%[text] LOS calculation using 1st order ROE mapping to RTN frame
%[text] ## inputs
%[text] `roe`: relative orbtial elements (ROE), nx6 matrix
%[text] $\delta \alpha =(a_d-a)/a \\\n\delta \lambda = (u_d - u ) + (\Omega_d - \Omega) \cos{i} \\\n\delta e_x = e_{xd} - e_{x} \\\n\delta e_y = e_{yd} - e_{y} \\\n\delta i_x = i_d - i \\\n\delta i_y = (\Omega_d - \Omega) \sin{i}\n$
%[text] `chiefOE`: absolute orbital elements of chief, nx6 matrix
%[text] $a , e, i, \Omega, w, f({\rm or} \~M)$ (km or m, -, rad, rad, rad, rad)
%[text] `flag`: 1 = true anomaly, 0:= mean anomaly , scalar
%[text] `GE`: gravitational constant of the Earth, m or km (unit must be unified with the position and velocity)
%[text] ## outputs
%[text] `azi`: R-T平面のR軸からの角度, nx1 vector, rad
%[text] `ele`: R-T平面からN方向への角度, nx1 vector, rad
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
