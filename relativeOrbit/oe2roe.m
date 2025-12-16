%[text] # calculating relative orbital elements (ROE) from absolute orbital elements
%[text] `chief`: absolute orbital elements of chief, nx6 matrix
%[text] $a , e, i, \\Omega, w, f({\\rm or} \~M)$ (km, -, rad, rad, rad, rad)
%[text] `deputy`: absolute orbital elements of deputy, nx6 matrix
%[text] $a , e, i, \\Omega, w, f({\\rm or} \~M)$ (km, -, rad, rad, rad, rad)
%[text] `roe`: relative orbtial elements (ROE), nx6 matrix
%[text] $\\delta \\alpha =(a\_d-a)/a \\\\\n\\delta \\lambda = (u\_d - u ) + (\\Omega\_d - \\Omega) \\cos{i} \\\\\n\\delta e\_x = e\_{xd} - e\_{x} \\\\\n\\delta e\_y = e\_{yd} - e\_{y} \\\\\n\\delta i\_x = i\_d - i \\\\\n\\delta i\_y = (\\Omega\_d - \\Omega) \\sin{i}\n$
%[text] `flag`: 1 = true anomaly, 0:= mean anomaly , scalar
%[text] ## note
%[text] NA
%[text] ## references
%[text] NA
%[text] ## revisions
%[text] 20211027  y.yoshimura, y.yoshimula@gmail.com
%[text] See also roe2DeputyOE.
function roe = oe2roe(chief, deputy, anomalyFlag)
% arguments
%     chief (:,6) {mustBeNumeric}
%     deputy (:,6) {mustBeNumeric}
%     anomalyFlag (1,1) {mustBeMember(anomalyFlag, [0, 1])}
% end

% range: [0, 2pi] for RAAN, w, f/M
chief(:,4:6) = mod(chief(:,4:6), 2*pi);
deputy(:,4:6) = mod(deputy(:,4:6), 2*pi);

if anomalyFlag == 1
    mC = meanAnomaly(chief(:,2), chief(:,6));
    mD = meanAnomaly(deputy(:,2), deputy(:,6));
else
    mC = chief(:,6);
    mD = deputy(:,6);
end

%[text] ## chief
aC = chief(:,1); % semi-major axis, km
uC = chief(:,5) + mC; % w + M, mean argument of latitude, rad
uC = mod(uC, 2 * pi);
exC = chief(:,2) .* cos(chief(:,5)); % e * cos(w), eccentric vector
eyC = chief(:,2) .* sin(chief(:,5)); % e * sin(w), eccentric vector
incC = chief(:,3); % inclination, rad
raanC = chief(:,4); % right ascension of the ascending node

%[text] ## deputy
aD = deputy(:,1); % semi-major axis, km
uD = deputy(:,5) + mD; % w + M, mean argument of latitude, rad
uD = mod(uD, 2 * pi);
exD = deputy(:,2) .* cos(deputy(:,5)); % e * cos(w), eccentric vector
eyD = deputy(:,2) .* sin(deputy(:,5)); % e * sin(w), eccentric vector
incD = deputy(:,3); % inclination, rad
raanD = deputy(:,4); % right ascension of the ascending node

%[text] ## ROEs
%[text] $\\delta \\alpha =(a\_d-a)/a \\\\\n\\delta \\lambda = (u\_d - u ) + (\\Omega\_d - \\Omega) \\cos{i} \\\\\n\\delta e\_x = e\_{xd} - e\_{x} \\\\\n\\delta e\_y = e\_{yd} - e\_{y} \\\\\n\\delta i\_x = i\_d - i \\\\\n\\delta i\_y = (\\Omega\_d - \\Omega) \\sin{i}\n$
deltaA = (aD - aC) ./ aC;
deltaLambda = uD - uC + (raanD - raanC) .* cos(incC);
deltaEx = exD - exC;
deltaEy = eyD - eyC;
deltaIx = incD - incC;
deltaIy = (raanD - raanC) .* sin(incC);

deltaLambda = wrapPi(deltaLambda); % wrap angle ( [-\pi, \pi) に変換 )

roe = [deltaA, deltaLambda, deltaEx, deltaEy, deltaIx, deltaIy];
end

%[appendix]{"version":"1.0"}
%---
