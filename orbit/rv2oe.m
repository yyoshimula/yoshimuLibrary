%[text] # position ${\\bf r}\n$ and velocity $\\bf v$ to orbital elements
%[text] `a`: semi-major axis, km
%[text] `e`: eccentricity
%[text] `inc`: inclination, rad
%[text] `raan`: longitude of the ascending node, rad
%[text] `w`: argument of perigee, rad
%[text] `nu`: ture anomaly, rad
%[text] `M`: mean anomaly at EPOCH, rad
%[text] `mu`: Earth's gravitational constant, scalar
%[text] ## note
%[text] NA
%[text] ## references 
%[text] Vallado, D.A., & Wayne D. McClain. Fundamentals of Astrodynamics and Applications. 4th edition, Springer Science & Business Media, 2001. pp114.
%[text] ## revisions
%[text] 20221110  y.yoshimura, y.yoshimula@gmail.com, y.yoshimura.a64@m.kyushu-u.ac.jp
%[text] See also oe2rv.
function oe = rv2oe(r, v, mu)
% arguments
%     r (:,3) {mustBeNumeric}
%     v (:,3) {mustBeNumeric}
%     mu
% end

nt_ = size(r, 1); % data length

%[text] ### orbital angular momentum, n vector and eccentricity vector
% orbital angular momentum
h = cross(r, v, 2); % 1x3 vector

nVec = cross([zeros(nt_,2), ones(nt_,1)], h, 2);

% eccentricity vector
eVec = (vecnorm(v,2,2).^2 - mu ./ vecnorm(r,2,2)) .* r - dot(r, v, 2) .* v;
eVec = eVec ./ mu;

% eccentricity
e = vecnorm(eVec,2,2);

% specfic mechanical energy, t_x1 vector
xi = vecnorm(v,2,2).^2 ./ 2 - mu ./ vecnorm(r,2,2);

% semi-major axis
paraInd = abs(xi) <= eps; % index when e = 1
a = (abs(xi) > eps) .* (-mu / 2 ./ xi);
a(paraInd) = Inf; % parabolic orbits

% semi-latus rectum
p = vecnorm(h, 2, 2).^2 ./ mu;

% inclination
inc = acos(h(:,3) ./ vecnorm(h,2,2));

%[text] ## type of orbits
%[text] 1 = elliptical and inclined
%[text] 2 = circular equatorial
%[text] 3 = circular inclined
%[text] 4 = elliptical, parabolic, hyperbolic equatorial

typeOrbit = ones(nt_,1); % elliptical and inclined
for i = 1:nt_
    if ( norm(eVec(i,:)) < eps)
        % ----------------  circular equatorial ---------------
        if  (inc(i)<eps) || (abs(inc(i)-pi) < eps)
            typeOrbit(i) = 2;
        else
            % --------------  circular inclined ---------------
            typeOrbit(i) = 3;
        end
    else
        % - elliptical, parabolic, hyperbolic equatorial --
        if  (inc(i) < eps) || (abs(inc(i) - pi) < eps)
            typeOrbit(i) = 4;
        end
    end
end
%[text] ### longitude of ascending node
raan = acos(nVec(:,1) ./ vecnorm(nVec,2,2));
raan = (nVec(:,2) >= 0.0) .* raan + (nVec(:,2) < 0.0) .* (2 * pi - raan);

%[text] ### true anomaly
nu = dot(r, eVec, 2) ./ vecnorm(r,2,2) ./ e;
nu = acos(nu);
nu = (dot(r,v,2) >= 0.0) .* nu + (dot(r,v,2) < 0.0) .* (2 * pi - nu);

%[text] ### argument of latitude for circular inlined and elliptical inclined orbits
u = dot(r, nVec, 2) ./ vecnorm(r,2,2) ./ vecnorm(nVec,2,2);
u = acos(u);
u = (r(:,3) >= 0.0) .* u + (r(:,3) < 0.0) .* (2 * pi - u);

%[text] ### longitude of perigee
for i = 1:nt_
    if typeOrbit(i) == 1
        % for elliptical inclined orbit
        w = dot(nVec, eVec, 2) ./ vecnorm(nVec,2,2) ./ e;
        w = acos(w);
        w = (eVec(:,3) >= 0.0) .* w + (eVec(:,3) < 0.0) .* (2 * pi - w);
    elseif typeOrbit(i) == 4 % ellptical equatorial
        wEE = (vecnorm(eVec,2,2) > eps) .* acos(eVec(:,1) ./ vecnorm(eVec,2,2));
        wEE = (eVec(:,2) < 0.0) .* (2*pi - wEE) + (eVec(:,2) >= 0.0) .* wEE;
        wEE = (inc >= pi / 2) .* (2*pi - wEE) + (inc < pi / 2) .* wEE;

        w = wEE;
    else

    end
end

%[text] ### ture longitude for circular equatorial
tmp = (vecnorm(r,2,2) > eps) .* r(:,1) ./ vecnorm(r,2,2);
trueLon = acos(tmp);
trueLon = (r(:,2) >= 0.0) .* trueLon + (r(:,2) < 0.0) .* (2*pi - trueLon);
trueLon = (inc < pi / 2) .* trueLon + (inc >= pi / 2) .* (2*pi - trueLon);

%[text] ## mean anomaly
M = meanAnomaly(e, nu);

%[text] ## orbital elements
oe = [a, e, inc, raan, w, nu];%, u, M, trueLon];

end

%[appendix]{"version":"1.0"}
%---
