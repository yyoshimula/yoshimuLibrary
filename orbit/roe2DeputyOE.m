%[text] # calculating deputy's absolute orbital elements using ROE and chief's absolute orbital elements
%[text] `roe`: relative orbtial elements (ROE), nx6 matrix
%[text] $\\delta \\alpha =(a\_d-a)/a \\\\\n\\delta \\lambda = (u\_d - u ) + (\\Omega\_d - \\Omega) \\cos{i} \\\\\n\\delta e\_x = e\_{xd} - e\_{x} \\\\\n\\delta e\_y = e\_{yd} - e\_{y} \\\\\n\\delta i\_x = i\_d - i \\\\\n\\delta i\_y = (\\Omega\_d - \\Omega) \\sin{i}\n$
%[text] `chiefOE`: absolute orbital elements of chief, nx6 matrix
%[text] $a , e, i, \\Omega, w, f({\\rm or} \~M)$ (km, -, rad, rad, rad, rad)
%[text] `flag`: 1 = true anomaly, 0:= mean anomaly , scalar
%[text] `deputyOE`: absolute orbital elements of deputy, nx6 matrix
%[text] $a , e, i, \\Omega, w, f({\\rm or} \~M)$ (km, -, rad, rad, rad, rad)
%[text] ## note
%[text] NA
%[text] ## references 
%[text] Di Mauro, G. 2019 Minimum-Fuel Control Strategy for Spacecraft Formation Reconfiguration via Finite-Time Maneuvers Journal of Guidance, Control, and Dynamics.
%[text] ## revisions
%[text] 20211027  y.yoshimura, y.yoshimula@gmail.com
%[text] See also oe2roe.
function deputyOE = roe2DeputyOE(roe, chiefOE, flag)
% arguments
%     roe (:,6) {mustBeNumeric}
%     chiefOE (:,6) {mustBeNumeric}
%     flag (1,1) {mustBeMember(flag, [0, 1])}
% end

%[text] ## relative orbital elements
dA = roe(:,1);
dLambda = roe(:,2);
dEx = roe(:,3);
dEy = roe(:,4);
dIx = roe(:,5);
dIy = roe(:,6);
%[text] ## chief (absolute) orbital elements
aC = chiefOE(:,1); % km
eC = chiefOE(:,2);
iC = chiefOE(:,3); % rad
raanC = chiefOE(:,4); % rad
raanC = mod(raanC, 2*pi);
wC = chiefOE(:,5); % rad
wC = mod(wC, 2*pi);
if flag == 1
    mC = meanAnomaly(eC, chiefOE(:,6)); % rad
else
    mC = chiefOE(:,6); % rad
end
mC = mod(mC, 2*pi);
%[text] ## deputy (absolute) orbital elements
%[text] $a\_{d} = a\\delta a + a \\\\\nu\_{d} = u + \\delta\\lambda - (\\Omega\_{d} - \\Omega)\\cos{i} \\\\\ne\_{d} = \\sqrt{(e\_{x}+\\delta e\_{x})^{2} + (e\_{y}+ \\delta e\_{y})^{2}} \\\\\ni\_{d} = i\_{c} + \\delta i\_{x} \\\\\n\\omega\_{d} = \\arctan{\\left(\\frac{e\_{y}+ \\delta e\_{y}}{e\_{x}+ \\delta e}\\right)} \\\\\n\\Omega\_{d} = \\Omega + \\frac{\\delta i\_{y}}{\\sin{i}}$
a = aC * dA + aC; % km
e = sqrt((dEx + eC .* cos(wC)).^2 + (dEy + eC .* sin(wC)).^2);

inc = iC + dIx;
raan = raanC + dIy ./ sin(iC);
raan = mod(raan, 2*pi);
w = atan2(dEy + eC .* sin(wC), dEx + eC .* cos(wC));
w = mod(w, 2*pi);
m = dLambda + wC + mC - (raan - raanC) .* cos(iC) - w;
m = mod(m , 2*pi);

if flag == 1 % true anomaly
    [f, ~] = trueAnomaly(a, e, m);
    deputyOE = [a, e, inc, raan, w, f];
else
    deputyOE = [a, e, inc, raan, w, m];
end


end

%[appendix]{"version":"1.0"}
%---
