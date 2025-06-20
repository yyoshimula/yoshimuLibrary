%[text] # Double-r method of orbit determination
%[text] `t:` time at observations, 3x1 vector
%[text] `azi, ele:` azimuth and elevation angles of the object at topocentric equatorial frame, 3x2 matrix
%[text] `rObs:` observer position vector at inertial frame, 3x3 matrix
%[text] `rRange`: initial guess range for the 8th order equations
%[text] `mu`: Earth's gravitational constant, scalar
%[text] ## note
%[text] NA
%[text] ## references 
%[text] Vallado, D.A., & Wayne D. McClain. Fundamentals of Astrodynamics and Applications. Springer Science & Business Media, 2001, pp443-445.
%[text] ## revisions
%[text] 20221202  y.yoshimura, y.yoshimula@gmail.com, y.yoshimura.a64@m.kyushu-u.ac.jp
%[text] See also gibbs.
function [r2, v2] = doubleR(t, aziele, rObs, mu, RE)

TOL = 1e-8 * RE; % tolerance for 8th order equation
r1Old = 9999.99;
r2Old = 9999.99;
%[text] ## time intervals
t1 = t(1);
t2 = t(2);
t3 = t(3);

tau1 = day2s(t1 - t2); % s
tau3 = day2s(t3 - t2);
%[text] ## initial guess
r(1) = 2.0 * RE; % km
r(2) = 2.01 * RE; % km
%[text] ## unit vectors
[L(:,1), L(:,2), L(:,3)] = sph2cart(aziele(:,1), aziele(:,2), 1); % line-of-sight direction, 3x3 matrix

% para
c = 2 .* dot(L, rObs, 2);
c(3,:) = []; % 2x1

%[text] ## iteration
while (abs(r(1) - r1Old) > TOL || abs(r(2) - r2Old) > TOL)

%[text] ### partial derivative w.r.t. $r\_1$
    [F1, F2, ~, ~, ~] = calcF(tau1, tau3, c, rObs, r, L, mu);
    %     Q = sqrt(F1^2 + F2^2);

    dltR1 = 0.005 * r(1);
    [dF1, dF2, ~, ~, ~] = calcF(tau1, tau3, c, rObs, [r(1)+dltR1, r(2)], L, mu);

    dF1dr1 = (dF1 - F1) / dltR1;
    dF2dr1 = (dF2 - F2) / dltR1;

%[text] ### partial derivative w.r.t. $r\_2$
    dltR2 = 0.005 * r(2);
    [dF1, dF2, ~, ~, ~] = calcF(tau1, tau3, c, rObs, [r(1), r(2)+dltR2], L, mu);

    dF1dr2 = (dF1 - F1) / dltR2;
    dF2dr2 = (dF2 - F2) / dltR2;

%[text] ### evaluate
    dlt = dF1dr1 * dF2dr2 - dF2dr1 * dF1dr2;

    dlt1 = dF2dr2 * F1 - dF1dr2 * F2;
    dlt2 = dF1dr1 * F2 - dF2dr1 * F1;

    dltR1 = - dlt1 / dlt;
    dltR2 = - dlt2 / dlt;

    r1Old = r(1) + dltR1;
    r2Old = r(2) + dltR2;

    % update
    r(1) = r(1) + dlt1;
    r(2) = r(2) + dlt2;    

end
%[text] ## iteration end
%[text] ### outpus: velocity vector
[~, ~, rVec, a, dltE] = calcF(tau1, tau3, c, rObs, r, L, mu);

f = 1 - a / r(2) * (1 - cos(dltE));
g = tau3 - sqrt(a^3 / mu) * (dltE - sin(dltE));

r2 = rVec(2,:);
v2 = (rVec(3,:) - f .* rVec(2,:)) ./ g;

end

%%
%[text] ## (local function) calculating F for partial derivatives
function [F1, F2, rVec, a, dltEout] = calcF(tau1, tau3, c, rObs, r, L, mu)
r = r(:);
c = c(:);
rho = -c + sqrt(c.^2 - 4 * (vecnorm(rObs(1:2,:),2,2).^2 - r(1:2,:).^2));
rho = rho ./ 2; % 2x1

rVec = rho .* L(1:2,:) + rObs(1:2,:); % 2x3 matrix

W = cross(rVec(1,:), rVec(2,:)); % 1x3 vector
W = W ./ norm(rVec(1,:)) ./ norm(rVec(2,:));

rho(3) = (-rObs(3,:) * W') / (L(3,:) * W');
rVec(3,:) = rho(3) .* L(3,:) + rObs(3,:);
r(3) = norm(rVec(3,:));

cosDltV = zeros(3,2);
sinDltV = zeros(3,2);
for j = 2:3
    for k = 1:2
        cosDltV(j,k) = (rVec(j,:) * rVec(k,:)') / norm(rVec(j,:)) / norm(rVec(k,:));
        sinDltV(j,k) = sqrt(1 - cosDltV(j,k)^2);
    end
end

dltV(2,1) = atan2(sinDltV(2,1), cosDltV(2,1));
dltV(3,1) = atan2(sinDltV(3,1), cosDltV(3,1));
if dltV(3,1) > pi
    c(1) = r(2) * sinDltV(3,2) / r(1) / sinDltV(3,1);
    c(3) = r(2) * sinDltV(2,1) / r(3) / sinDltV(3,1);
    p = (c(1) * r(1) + c(3) * r(3) - r(2)) / (c(1) + c(3) - 1);
else
    c(1) = r(1) * sinDltV(3,1) / r(2) / sinDltV(3,2);
    c(3) = r(1) * sinDltV(2,1) / r(3) / sinDltV(3,2);
    p = (c(3) * r(3) - c(1) * r(2) + r(1)) / (-c(1) + c(3) + 1);
end

ecosV = p ./ r - 1; % 3x1 vector

if dltV(2,1) <= deg2rad(179) || deg2rad(181) <= dltV(2,1)
    esinV(2) = (-cosDltV(2,1) * ecosV(2) + ecosV(1)) / sinDltV(2,1);
else
    esinV(2) = (cosDltV(3,2) * ecosV(2) - ecosV(3)) / sinDltV(3,2); %% textbook has error, 32 is correct, not 31
end

e = sqrt(ecosV(2)^2 + esinV(2)^2);

a = p / (1 - e^2);
n = sqrt(mu / a^3);

S = r(2) / p * sqrt(1 - e^2) * esinV(2);
C = r(2) / p * (e^2 + ecosV(2));

sinDltE(3,2) = r(3) / sqrt(a * p) * sinDltV(3,2) - r(3) / p * (1 - cosDltV(3,2)) * S;
cosDltE(3,2) = 1 - r(2) * r(3) / (a * p) * (1 - cosDltV(3,2));

sinDltE(2,1) = r(1) / sqrt(a * p) * sinDltV(2,1) + r(1) / p * (1 - cosDltV(2,1)) * S;
cosDltE(2,1) = 1 - r(2) * r(1) / (a * p) * (1 - cosDltV(2,1));

dltE(3,2) = atan2(sinDltE(3,2), cosDltE(3,2));
dltE(2,1) = atan2(sinDltE(2,1), cosDltE(2,1));

dltM(3,2) = dltE(3,2) + 2 * S * sin(dltE(3,2)/2)^2 - C * sin(dltE(3,2));
dltM(1,2) = -dltE(2,1) + 2 * S * sin(dltE(2,1)/2)^2 + C * sin(dltE(2,1));

%[text] ### outputs
F1 = tau1 - dltM(1,2) / n;
F2 = tau3 - dltM(3,2) / n;

dltEout = dltE(3,2);

end

%[appendix]{"version":"1.0"}
%---
