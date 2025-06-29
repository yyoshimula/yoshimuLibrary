%[text] # mapping relative orbital elements (ROEs) to relative position and velocity at RTN frame
%[text] ## inputs
%[text] `roe`: relative orbtial elements (ROE), nx6 matrix
%[text] $\\delta \\alpha =(a\_d-a)/a \\\\\n\\delta \\lambda = (u\_d - u ) + (\\Omega\_d - \\Omega) \\cos{i} \\\\\n\\delta e\_x = e\_{xd} - e\_{x} \\\\\n\\delta e\_y = e\_{yd} - e\_{y} \\\\\n\\delta i\_x = i\_d - i \\\\\n\\delta i\_y = (\\Omega\_d - \\Omega) \\sin{i}\n$
%[text] `chiefOE`: absolute orbital elements of chief, nx6 matrix
%[text] $a , e, i, \\Omega, w, f({\\rm or} \~M)$ (m or km, -, rad, rad, rad, rad)
%[text] `flag`: 1 = true anomaly, 0:= mean anomaly , scalar
%[text] `GE`: gravitational constant of the Earth, m or km (unit must be unified with the position and velocity)
%[text] ## outputs
%[text] xRTN: deputy position @RTN frame, m or km
%[text] vRTN: deputy velocity @RTN frame, m or km/s
%[text] ## note
%[text] 重力定数と位置・速度は単位を合わせること
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20211027  y.yoshimura, y.yoshimula@gmail.com
%[text] See also .
function [xRTN, vRTN] = roe2rtn(roe, chiefOE, flag, GE)
% arguments
%     roe (:,6) {mustBeNumeric}
%     chiefOE (:,6) {mustBeNumeric}
%     flag (1,1) {mustBeMember(flag, [0, 1])}
%     GE (1,1) {mustBeNumeric}
% end

roe(:,2) = atan2(sin(roe(:,2)), cos(roe(:,2)));

chiefOE(:,6) = mod(chiefOE(:,6), 2*pi);
aC = chiefOE(:,1); % km, semi-major axis of chief
nC = sqrt(GE ./ aC.^3); % mean motion of chief

if flag == 1
   mC = meanAnomaly(chiefOE(:,2), chiefOE(:,6));
else
   mC = chiefOE(:,6);
end
mC = mod(mC, 2*pi);
uC = chiefOE(:,5) + mC; % mean argument of latitude, rad
uC = mod(uC, 2*pi);
%[text] ## mapping
x = roe(:,1) - cos(uC) .* roe(:,3) - sin(uC) .* roe(:,4);
y = roe(:,2) + 2 .* sin(uC) .* roe(:,3) - 2 .* cos(uC) .* roe(:,4);

vx = nC .* sin(uC) .* roe(:,3) - nC .* cos(uC) .* roe(:,4);
vy = -3/2 .* nC .* roe(:,1) + 2 .* nC .* cos(uC) .* roe(:,3) + 2 * nC .* sin(uC) .* roe(:,4);

z = sin(uC) .* roe(:,5) - cos(uC) .* roe(:,6);
vz = nC .* cos(uC) .* roe(:,5) + nC .* sin(uC) .* roe(:,6);

xRTN = aC .* [x y z];
vRTN = aC .* [vx vy vz];

end

%[appendix]{"version":"1.0"}
%---
