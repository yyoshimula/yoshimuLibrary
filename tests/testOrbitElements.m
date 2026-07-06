classdef testOrbitElements < matlab.unittest.TestCase
    % Unit tests for orbital-element utilities in orbit/:
    %   keplerEq, oe2rv, rv2oe, coe2mee, mee2coe, trueAnomaly,
    %   meanAnomaly, geodetic2Geocentric, geocentric2Geodetic, orbitConst.
    %
    % Units follow the library conventions: km, km/s, rad.
    %
    % Note: rv2oe is intentionally NOT tested with circular orbits
    % (typeOrbit 2 and 3): for those inputs the internal variable `w`
    % is never assigned (orbit/rv2oe.m, loop at lines 91-106) and the
    % function errors at line 118. Only non-circular inclined orbits
    % are exercised here.

    methods (TestClassSetup)
        function addLibraryToPath(testCase) %#ok<MANU>
            root = fileparts(fileparts(mfilename('fullpath')));
            % genpath includes hidden folders such as .claude/worktrees,
            % whose stale library copies would shadow the real files
            p = strsplit(genpath(root), pathsep);
            p = p(~cellfun(@isempty, p) & ~contains(p, [filesep '.']));
            addpath(strjoin(p, pathsep));
        end
    end

    methods (Test)
        %% keplerEq -------------------------------------------------------
        function testKeplerEqResidualSweep(testCase)
            % Solution E of Kepler's equation must satisfy
            % E - e*sin(E) - M = 0 to well below 1e-12 (Newton iteration
            % terminates on a step <= 1e-8; quadratic convergence leaves
            % a residual of order eps).
            for e = 0:0.1:0.9
                for M = 0:pi/12:(2*pi - 1e-9)
                    E = keplerEq(M, e);
                    raw = E - e*sin(E) - M;
                    % keplerEq wraps E into [0, 2*pi]; remove any 2*pi
                    % offset before measuring the residual
                    res = abs(raw - 2*pi*round(raw/(2*pi)));
                    testCase.verifyEqual(res, 0, 'AbsTol', 1e-12, ...
                        sprintf('Kepler residual too large for e=%.1f, M=%.6f', e, M));
                end
            end
        end

        function testKeplerEqZeroEccentricity(testCase)
            % For e = 0 the eccentric anomaly equals the mean anomaly.
            for M = [0, 0.5, pi/2, pi, 4.0, 6.0]
                E = keplerEq(M, 0);
                testCase.verifyEqual(E, M, 'AbsTol', 1e-12);
            end
        end

        %% oe2rv ----------------------------------------------------------
        function testOe2rvCircularEquatorial(testCase)
            % A circular equatorial orbit with nu = 0 must give
            % r = [a 0 0] and v = [0 sqrt(mu/a) 0] (pure algebra).
            const = orbitConst();
            a = 7000; % km
            [r, v] = oe2rv([a, 0, 0, 0, 0, 0], 1, const.GE);
            vCirc = sqrt(const.GE / a); % km/s
            testCase.verifyEqual(r, [a, 0, 0], 'AbsTol', 1e-9);
            testCase.verifyEqual(v, [0, vCirc, 0], 'AbsTol', 1e-12);
        end

        function testOe2rvMeanAnomalyFlagConsistency(testCase)
            % flag = 0 (mean anomaly) must be equivalent to converting
            % M -> nu with trueAnomaly and calling oe2rv with flag = 1.
            const = orbitConst();
            a = 8000; e = 0.2; inc = 1.0; raan = 0.7; w = 2.5; M = 1.3;
            [rM, vM] = oe2rv([a, e, inc, raan, w, M], 0, const.GE);
            nu = trueAnomaly(a, e, M);
            [rNu, vNu] = oe2rv([a, e, inc, raan, w, nu], 1, const.GE);
            testCase.verifyEqual(rM, rNu, 'AbsTol', 1e-6);  % km
            testCase.verifyEqual(vM, vNu, 'AbsTol', 1e-9);  % km/s
        end

        %% rv2oe <-> oe2rv round trip -------------------------------------
        function testRv2oeOe2rvRoundTripEllipticInclined(testCase)
            % oe -> (r, v) -> oe must reproduce the classical elements
            % [a, e, inc, raan, w, nu] for generic non-circular,
            % non-equatorial orbits.
            const = orbitConst();
            oeSet = [8000,  0.1, 0.5,            1.0, 2.0, 3.0;   % LEO-ish
                     26560, 0.7, deg2rad(63.4),  4.0, 5.0, 0.5];  % Molniya-like
            for k = 1:size(oeSet, 1)
                oe0 = oeSet(k, :);
                [r, v] = oe2rv(oe0, 1, const.GE);
                oe1 = rv2oe(r, v, const.GE);
                testCase.verifyEqual(oe1(1), oe0(1), 'RelTol', 1e-9, ...
                    'semi-major axis a not recovered');
                testCase.verifyEqual(oe1(2), oe0(2), 'AbsTol', 1e-12, ...
                    'eccentricity e not recovered');
                testCase.verifyEqual(oe1(3:6), oe0(3:6), 'AbsTol', 1e-9, ...
                    'angles [inc, raan, w, nu] not recovered');
            end
        end

        function testRv2oeCircularEquatorialHandBuilt(testCase)
            % Regression: circular orbits used to error with
            % "Unrecognized function or variable 'w'" because the
            % typeOrbit 2/3 branch never assigned w.
            % For a hand-built circular equatorial state the convention is
            % w = 0, raan = 0, nu = true longitude (angle of r from +x).
            const = orbitConst();
            a = 7000; % km
            vc = sqrt(const.GE / a); % km/s
            % state at the +x axis (the reported failing input)
            oe = rv2oe([a, 0, 0], [0, vc, 0], const.GE);
            testCase.verifyEqual(oe, [a, 0, 0, 0, 0, 0], 'AbsTol', 1e-9);
            % state at true longitude 2.0 rad
            th = 2.0;
            oe = rv2oe(a*[cos(th), sin(th), 0], vc*[-sin(th), cos(th), 0], ...
                const.GE);
            testCase.verifyEqual(oe(1), a, 'RelTol', 1e-12);
            testCase.verifyEqual(oe(2), 0, 'AbsTol', 1e-12);
            testCase.verifyEqual(oe(3:5), [0, 0, 0], 'AbsTol', 1e-12);
            testCase.verifyEqual(oe(6), th, 'AbsTol', 1e-9);
        end

        function testRv2oeCircularPolarHandBuilt(testCase)
            % Circular polar orbit crossing the north pole: raan = 0,
            % w = 0, nu = argument of latitude = pi/2 (exact geometry:
            % r along +z, node line along +x).
            const = orbitConst();
            a = 7000; % km
            vc = sqrt(const.GE / a); % km/s
            oe = rv2oe([0, 0, a], [-vc, 0, 0], const.GE);
            testCase.verifyEqual(oe(1), a, 'RelTol', 1e-12);
            testCase.verifyEqual(oe(2), 0, 'AbsTol', 1e-12);
            testCase.verifyEqual(oe(3), pi/2, 'AbsTol', 1e-12);
            testCase.verifyEqual(oe(4), 0, 'AbsTol', 1e-12);
            testCase.verifyEqual(oe(5), 0, 'AbsTol', 1e-12);
            testCase.verifyEqual(oe(6), pi/2, 'AbsTol', 1e-9);
        end

        function testRv2oeOe2rvRoundTripCircularInclined(testCase)
            % oe -> (r, v) -> oe for a circular inclined orbit: rv2oe
            % returns w = 0 and nu = argument of latitude. With w = 0 in
            % the input, the angle from the ascending node equals nu, so
            % [a, 0, inc, raan, 0, nu] must be recovered exactly.
            const = orbitConst();
            oe0 = [7000, 0, deg2rad(51.6), 1.0, 0, 2.5];
            [r, v] = oe2rv(oe0, 1, const.GE);
            oe1 = rv2oe(r, v, const.GE);
            testCase.verifyEqual(oe1(1), oe0(1), 'RelTol', 1e-9);
            testCase.verifyEqual(oe1(2), 0, 'AbsTol', 1e-12);
            testCase.verifyEqual(oe1(3:6), oe0(3:6), 'AbsTol', 1e-9);
        end

        function testRv2oeOe2rvRoundTripEllipticEquatorial(testCase)
            % Elliptical equatorial orbit (typeOrbit 4): raan used to come
            % back NaN (zero node vector). Convention: raan = 0 and w is
            % the true longitude of perigee.
            const = orbitConst();
            oe0 = [8000, 0.1, 0, 0, 1.2, 0.7];
            [r, v] = oe2rv(oe0, 1, const.GE);
            oe1 = rv2oe(r, v, const.GE);
            testCase.verifyEqual(oe1(1), oe0(1), 'RelTol', 1e-9);
            testCase.verifyEqual(oe1(2), oe0(2), 'AbsTol', 1e-12);
            testCase.verifyEqual(oe1(3:6), oe0(3:6), 'AbsTol', 1e-9);
        end

        function testRv2oeVectorizedMixedTypes(testCase)
            % One call with four stacked states, one of each orbit type,
            % must recover every row (the old per-type loop overwrote w
            % for all rows on each iteration).
            const = orbitConst();
            oeSet = [8000, 0.1, 0.5,            1.0, 2.0, 3.0;  % elliptical inclined
                     7000, 0,   0,              0,   0,   2.0;  % circular equatorial
                     7000, 0,   deg2rad(51.6),  1.0, 0,   2.5;  % circular inclined
                     8000, 0.1, 0,              0,   1.2, 0.7]; % elliptical equatorial
            n = size(oeSet, 1);
            r = zeros(n, 3); v = zeros(n, 3);
            for k = 1:n
                [r(k,:), v(k,:)] = oe2rv(oeSet(k,:), 1, const.GE);
            end
            oe1 = rv2oe(r, v, const.GE);
            testCase.verifyEqual(oe1(:,1), oeSet(:,1), 'RelTol', 1e-9, ...
                'semi-major axis a not recovered');
            testCase.verifyEqual(oe1(:,2), oeSet(:,2), 'AbsTol', 1e-12, ...
                'eccentricity e not recovered');
            testCase.verifyEqual(oe1(:,3:6), oeSet(:,3:6), 'AbsTol', 1e-9, ...
                'angles [inc, raan, w, nu] not recovered');
        end

        function testRv2oeVisVivaConsistency(testCase)
            % a returned by rv2oe must satisfy the vis-viva equation
            % xi = -mu/(2a) for an independent hand-built state.
            const = orbitConst();
            r = [6524.834, 6862.875, 6448.296];   % km (Vallado ex. 2-5 state)
            v = [4.901327, 5.533756, -1.976341];  % km/s
            oe = rv2oe(r, v, const.GE);
            xi = norm(v)^2/2 - const.GE/norm(r);
            aVisViva = -const.GE/(2*xi);
            testCase.verifyEqual(oe(1), aVisViva, 'RelTol', 1e-12);
            % specific angular momentum h = sqrt(mu*a*(1-e^2))
            hFromElems = sqrt(const.GE * oe(1) * (1 - oe(2)^2));
            testCase.verifyEqual(norm(cross(r, v)), hFromElems, 'RelTol', 1e-9);
        end

        %% coe2mee / mee2coe ----------------------------------------------
        function testCoe2MeeKnownValues(testCase)
            % Modified equinoctial elements from their definitions:
            % p = a(1-e^2), f = e*cos(w+raan), g = e*sin(w+raan),
            % h = tan(i/2)*cos(raan), k = tan(i/2)*sin(raan), L = raan+w+nu.
            oe.a = 8000; oe.e = 0.1; oe.inc = 0.5;
            oe.raan = 1.0; oe.w = 2.0; oe.nu = 3.0;
            mee = coe2mee(oe);
            testCase.verifyEqual(mee.p_, oe.a*(1 - oe.e^2),        'AbsTol', 1e-9);
            testCase.verifyEqual(mee.f_, oe.e*cos(oe.w + oe.raan), 'AbsTol', 1e-12);
            testCase.verifyEqual(mee.g_, oe.e*sin(oe.w + oe.raan), 'AbsTol', 1e-12);
            testCase.verifyEqual(mee.h_, tan(oe.inc/2)*cos(oe.raan), 'AbsTol', 1e-12);
            testCase.verifyEqual(mee.k_, tan(oe.inc/2)*sin(oe.raan), 'AbsTol', 1e-12);
            testCase.verifyEqual(mee.L_, mod(oe.raan + oe.w + oe.nu, 2*pi), 'AbsTol', 1e-12);
        end

        function testMeeCoeRoundTrip(testCase)
            % coe -> mee -> coe must be the identity for a generic
            % elliptic inclined orbit (angles already in [0, 2*pi)).
            oe0.a = 8000; oe0.e = 0.1; oe0.inc = 0.5;
            oe0.raan = 1.0; oe0.w = 2.0; oe0.nu = 3.0;
            oe1 = mee2coe(coe2mee(oe0));
            testCase.verifyEqual(oe1.a,    oe0.a,    'RelTol', 1e-12);
            testCase.verifyEqual(oe1.e,    oe0.e,    'AbsTol', 1e-12);
            testCase.verifyEqual(oe1.inc,  oe0.inc,  'AbsTol', 1e-12);
            testCase.verifyEqual(oe1.raan, oe0.raan, 'AbsTol', 1e-12);
            testCase.verifyEqual(oe1.w,    oe0.w,    'AbsTol', 1e-12);
            testCase.verifyEqual(oe1.nu,   oe0.nu,   'AbsTol', 1e-12);
        end

        %% trueAnomaly ----------------------------------------------------
        function testTrueAnomalyRoundTripAndRadius(testCase)
            % trueAnomaly(a, e, M) -> f must invert through
            % meanAnomaly(e, f) -> M, and the returned radius must equal
            % both a*(1 - e*cos(E)) and the conic equation p/(1+e*cos f).
            a = 8000; e = 0.3;
            for M = [0.1, 1.0, 2.0, pi, 4.5, 6.0]
                [f, r] = trueAnomaly(a, e, M);
                testCase.verifyEqual(meanAnomaly(e, f), M, 'AbsTol', 1e-9);
                E = keplerEq(M, e);
                testCase.verifyEqual(r, a*(1 - e*cos(E)), 'AbsTol', 1e-8);
                p = a*(1 - e^2);
                testCase.verifyEqual(r, p/(1 + e*cos(f)), 'AbsTol', 1e-6);
            end
        end

        function testTrueAnomalyCircular(testCase)
            % For e = 0: f = M and r = a.
            a = 7000; M = 2.2;
            [f, r] = trueAnomaly(a, 0, M);
            testCase.verifyEqual(f, M, 'AbsTol', 1e-12);
            testCase.verifyEqual(r, a, 'AbsTol', 1e-9);
        end

        %% geodetic <-> geocentric ----------------------------------------
        function testGeodeticGeocentricRoundTrip(testCase)
            % lat/lon/h -> ECEF -> lat/lon/h must be the identity.
            % geocentric2Geodetic iterates the latitude to 1e-8 rad with a
            % contraction factor ~e^2 (~0.0067), so the converged latitude
            % is good to well below 1e-9 rad; longitude is exact (atan2).
            const = orbitConst();
            latList = deg2rad([35, -60, 5, 80]);
            lonList = deg2rad([139, -100, 0, 170]);
            hList = [0.5, 2.3, 0, 10]; % km
            for k = 1:numel(latList)
                r = geodetic2Geocentric(latList(k), lonList(k), hList(k), ...
                    const.RE, const.fE);
                [lon2, lat2, h2] = geocentric2Geodetic(r(1), r(2), r(3), ...
                    const.RE, const.fE);
                testCase.verifyEqual(lat2, latList(k), 'AbsTol', 1e-9);
                testCase.verifyEqual(lon2, lonList(k), 'AbsTol', 1e-9);
                testCase.verifyEqual(h2,   hList(k),   'AbsTol', 1e-6); % km
            end
        end

        function testGeodeticEquatorFixedPoint(testCase)
            % On the equator at h = 0 the point is [RE 0 0] and the
            % inverse conversion returns (lon, lat, h) = (0, 0, 0).
            const = orbitConst();
            r = geodetic2Geocentric(0, 0, 0, const.RE, const.fE);
            testCase.verifyEqual(r, [const.RE, 0, 0], 'AbsTol', 1e-9);
            [lon, lat, h] = geocentric2Geodetic(const.RE, 0, 0, ...
                const.RE, const.fE);
            testCase.verifyEqual(lon, 0, 'AbsTol', 1e-12);
            testCase.verifyEqual(lat, 0, 'AbsTol', 1e-12);
            testCase.verifyEqual(h,   0, 'AbsTol', 1e-9);
        end

        function testGeodeticPoleForward(testCase)
            % At the north pole the forward conversion must give
            % [0, 0, b + h] with polar radius b = RE*(1 - f).
            % (The iterative inverse is not defined at the exact pole,
            % so only the forward direction is checked here.)
            const = orbitConst();
            h = 1.0; % km
            r = geodetic2Geocentric(pi/2, 0, h, const.RE, const.fE);
            b = const.RE*(1 - const.fE);
            testCase.verifyEqual(r(1), 0,     'AbsTol', 1e-9);
            testCase.verifyEqual(r(2), 0,     'AbsTol', 1e-9);
            testCase.verifyEqual(r(3), b + h, 'AbsTol', 1e-9);
        end

        %% orbitConst -----------------------------------------------------
        function testOrbitConstValues(testCase)
            % Literal constants (Vallado 4th ed. / WGS-84).
            const = orbitConst();
            testCase.verifyEqual(const.GE, 398600.4415, 'AbsTol', 0);
            testCase.verifyEqual(const.RE, 6378.137,    'AbsTol', 0);
            testCase.verifyEqual(const.J2000, 2451545.0, 'AbsTol', 0);
            testCase.verifyEqual(const.c, 299792458,    'AbsTol', 0);
            % derived values
            testCase.verifyEqual(const.fE, 1/298.257,   'RelTol', 1e-15);
            testCase.verifyEqual(const.GEm, const.GE*1e9, 'RelTol', 1e-15);
            testCase.verifyEqual(const.REm, const.RE*1e3, 'RelTol', 1e-15);
        end
    end
end
