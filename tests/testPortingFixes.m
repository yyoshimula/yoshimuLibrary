classdef testPortingFixes < matlab.unittest.TestCase
    % testPortingFixes  Regression tests for runtime bugs found while comparing the library
    % with its Julia port (YoshimuLibrary.jl, test/matlab/gen_*.m), fixed 2026-10-01.
    %
    %   - sclerp: failed for its documented 1x8 row inputs
    %   - roe2DeputyOE: failed for n x 6 inputs (matrix product aC * dA)
    %   - srpApproxCT2: lam = 2 / mCT.^2 (matrix right division) gave NaN for several facets
    %   - srpApproxCT / srpApproxCT2: called the nonexistent km2AU (the file is km2au.m)
    %   - jr1971: needed wrapToPi from the Mapping Toolbox (now the library's wrapPi)
    %
    % Found while syncing the Python port (yoshimulib), fixed 2026-10-05:
    %   - jr1971: the factor f was missing in the exponent between 90 and 100 km
    %     (almost constant density and a jump of about 6 at 100 km)
    %   - oe2roe / roe2DeputyOE: RAAN difference of a pair on both sides of 0/2pi
    %   - srpAS / srpASuni: the diffuse term was [0, 0, sum(cd1) ...] for every facet
    %     (valid only for a single facet whose normal is +z)
    %   - calcLocalFrame, calcNormalObj, srp* (torque), calcRelPosVelAtti: cross without
    %     dim worked along the columns for exactly 3 rows (3 facets / 3 time steps)
    %
    % Run headless with:
    %   matlab -batch "r = runtests('tests/testPortingFixes.m'); disp(table(r)); assertSuccess(r)"

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

    methods (Static, Access = private)
        function [dq1, dq2] = twoPoses()
            q1 = [0.1 -0.2 0.3 0.9]; q1 = q1 / norm(q1);
            q2 = [0.4 0.1 -0.2 0.8]; q2 = q2 / norm(q2);
            dq1 = pos2dq(1, 4, [1 -2 0.5], q1);
            dq2 = pos2dq(1, 4, [3 0.5 2], q2);
        end

        function [chief, deputy] = formation()
            chief = [7000, 1e-3, deg2rad(50), 0.7, 0.4, 1.1;
                     7100, 0.01, 1.2,         3.0, 2.0, 5.0;
                     6900, 0.02, 0.3,         5.0, 4.0, 2.0];
            deputy = chief + [ 0.1,   1e-4,  1e-4,  2e-4, -3e-3,  2e-3;
                              -0.2,  -2e-4,  2e-4, -1e-4,  1e-3, -1e-3;
                               0.05,  1e-4, -1e-4,  1e-4,  2e-3,  1e-3];
        end

        function sat = threeFacets(idx)
            % three triangular facets with different normals and optical properties
            normal = [0 0 1; 0.6 0 0.8; 0 -0.8 0.6];
            u = [1 0 0; 0 1 0; 1 0 0]; % in-plane edge directions
            sat.vertices = zeros(9, 3);
            for i = 1:3
                p0 = [i, 0, 0];
                sat.vertices(3*i-2:3*i, :) = [p0; p0 + u(i,:); p0 + cross(normal(i,:), u(i,:))];
            end
            sat.faces = [1 2 3; 4 5 6; 7 8 9];
            sat.normal = normal;
            sat.area = [1; 0.5; 2];
            sat.pos = [0 0 0; 1 0 0.3; 0 1 -0.2];
            sat.Cd = [0.3; 0.5; 0.1];
            sat.F0 = [0.5; 0.2; 0.7];
            sat.nu = [8; 40; 20];
            sat.nv = [40; 8; 20];
            if nargin > 0
                sat.faces = sat.faces(idx, :);
                for name = {'normal', 'area', 'pos', 'Cd', 'F0', 'nu', 'nv'}
                    sat.(name{1}) = sat.(name{1})(idx, :);
                end
            end
            [sat.uu, sat.uv, sat.qlb] = calcLocalFrame(sat);
        end

        function s = addRow(s)
            % one more time step: histories with 3 rows get a 4th row
            for name = fieldnames(s)'
                v = s.(name{1});
                if size(v, 1) == 3
                    s.(name{1}) = [v; v(end, :)];
                end
            end
        end
    end

    methods (Test)
        function sclerpRowInputs(testCase)
            [dq1, dq2] = testPortingFixes.twoPoses();
            t = (0:0.25:1)';
            dqt = sclerp(t, 4, dq1, dq2);
            testCase.verifySize(dqt, [numel(t), 8]);
            testCase.verifyEqual(dqt(1,:), dq1, 'AbsTol', 1e-12);
            testCase.verifyEqual(dqt(end,:), dq2, 'AbsTol', 1e-12);
            % scalar t gives the same rows; real parts stay unit quaternions
            for k = 1:numel(t)
                testCase.verifyEqual(sclerp(t(k), 4, dq1, dq2), dqt(k,:), 'AbsTol', 1e-14);
            end
            testCase.verifyEqual(vecnorm(dqt(:,1:4), 2, 2), ones(numel(t), 1), 'AbsTol', 1e-12);
        end

        function sclerpScalarFirst(testCase)
            [dq1, dq2] = testPortingFixes.twoPoses();
            s = @(dq) dq(:, [4 1 2 3 8 5 6 7]);          % scalar-last -> scalar-first
            dqt4 = sclerp([0; 0.4; 1], 4, dq1, dq2);
            dqt0 = sclerp([0; 0.4; 1], 0, s(dq1), s(dq2));
            testCase.verifyEqual(dqt0, s(dqt4), 'AbsTol', 1e-12);
        end

        function roe2DeputyOEMultiRow(testCase)
            [chief, deputy] = testPortingFixes.formation();
            for flag = [1, 0]
                roe = oe2roe(chief, deputy, flag);
                dep = roe2DeputyOE(roe, chief, flag);
                testCase.verifySize(dep, size(chief));
                for k = 1:size(chief, 1)
                    testCase.verifyEqual(dep(k,:), roe2DeputyOE(roe(k,:), chief(k,:), flag), ...
                        'AbsTol', 1e-12);
                end
                % round trip recovers the deputy elements (angles mod 2*pi)
                d = dep - deputy;
                d(:, 4:6) = mod(d(:, 4:6) + pi, 2*pi) - pi;
                testCase.verifyLessThan(max(abs(d(:, 1))), 1e-8);
                testCase.verifyLessThan(max(abs(d(:, 2:6)), [], 'all'), 1e-10);
            end
        end

        function roe2DeputyOEEquatorialRow(testCase)
            % an equatorial chief in one row must not change the other rows
            [chief, deputy] = testPortingFixes.formation();
            roe = oe2roe(chief, deputy, 0);
            chief2 = chief; chief2(2, 3) = 0;
            dep = roe2DeputyOE(roe, chief2, 0);
            testCase.verifyEqual(dep(2, 4), mod(chief2(2, 4), 2*pi), 'AbsTol', 1e-15);
            testCase.verifyEqual(dep([1 3], :), roe2DeputyOE(roe([1 3], :), chief([1 3], :), 0), ...
                'AbsTol', 1e-12);
        end

        function srpApproxCT2EqualsSumOfFacets(testCase)
            const = orbitConst;
            d = au2km(1.0, const) * 10^3; % m
            thetaN = deg2rad([0; 10; 30; 50; 75]);
            sat2.area = [1; 0.5; 2; 1; 0.3];
            sat2.F0 = [0.5; 0.4; 0.6; 0.5; 0.3];
            sat2.mCT = [0.1; 0.2; 0.3; 0.15; 0.25];
            srp2 = srpApproxCT2(sat2, thetaN, [0 0 1], d, const);
            testCase.verifyTrue(all(isfinite(srp2)));
            total = zeros(1, 3);
            for i = 1:numel(thetaN)
                s.area = sat2.area(i); s.F0 = sat2.F0(i); s.mCT = sat2.mCT(i);
                total = total + srpApproxCT(s, thetaN(i), [0 0 1], d, const);
            end
            testCase.verifyEqual(srp2, total, 'RelTol', 1e-12, 'AbsTol', 1e-18);
        end

        function jr1971WithoutMappingToolbox(testCase)
            % jr1971 now uses the library's wrapPi; densities stay finite and positive
            jd = gc2jd(2017, 1, 1, 0, 0, 0);
            for h = [120e3, 300e3, 800e3]
                out = jr1971(jd, deg2rad(20), deg2rad(60), h, 150, 100, 4);
                testCase.verifyTrue(isfinite(out.total_density) && out.total_density > 0);
            end
        end

        function jr1971Between90And100km(testCase)
            % barometric equation: the density scale height is about 5.6 km here and the
            % 90-100 km branch joins the 100-125 km branch at 100 km
            jd = gc2jd(2017, 1, 1, 0, 0, 0);
            cases = [20, 60, 150, 100, 4; -45, 200, 70, 75, 1; 80, -30, 220, 180, 7.3];
            for i = 1:size(cases, 1)
                c = cases(i, :);
                rho = @(hkm) getfield(jr1971(jd, deg2rad(c(1)), deg2rad(c(2)), hkm * 1e3, ...
                    c(3), c(4), c(5)), 'total_density'); %#ok<GFLD>
                H = -10 / log(rho(100) / rho(90)); % km
                testCase.verifyGreaterThan(H, 5);
                testCase.verifyLessThan(H, 6.5);
                testCase.verifyEqual(rho(100), rho(100 + 1e-6), 'RelTol', 1e-3);
            end
            % hydrostatic integration of the model's own T(z) and M(z) from 90 km
            % (Kp = 4, F10 = 150, F10a = 100, 20 deg N, 60 deg E, 2017-01-01 + 77.3 d).
            % The closed form cancels terms of about 3e6 in the exponent, so it is
            % accurate to about 1e-5 in double precision
            out = jr1971(jd + 77.3, deg2rad(20), deg2rad(60), 95e3, 150, 100, 4);
            testCase.verifyEqual(out.total_density, 1.609410e-06, 'RelTol', 1e-4);
        end

        function oe2roeRaanAcrossZero(testCase)
            % chief RAAN just below 2*pi, deputy RAAN just above 0 (and vice versa)
            for flag = [1, 0]
                for s = [1, -1]
                    chief = [7000, 1e-3, 0.9, 2*pi - s * 1e-4, 0.4, 1.1];
                    deputy = chief + [0.1, 1e-4, 1e-4, s * 3e-4, -3e-3, 2e-3];
                    deputy(4) = mod(deputy(4), 2*pi);
                    % the same pair rotated about the z-axis, away from 0/2*pi
                    shift = [0, 0, 0, 1, 0, 0];
                    roe = oe2roe(chief, deputy, flag);
                    testCase.verifyEqual(roe, oe2roe(chief - shift, deputy - shift, flag), ...
                        'AbsTol', 1e-12);
                    testCase.verifyEqual(roe(6), s * 3e-4 * sin(0.9), 'AbsTol', 1e-12);
                    % round trip (roe2DeputyOE wraps the deputy RAAN at 2*pi)
                    d = roe2DeputyOE(roe, chief, flag) - deputy;
                    d(4:6) = mod(d(4:6) + pi, 2*pi) - pi;
                    testCase.verifyLessThan(abs(d(1)), 1e-8);
                    testCase.verifyLessThan(max(abs(d(2:6))), 1e-10);
                end
            end
        end

        function srpASDiffuseTermPerFacet(testCase)
            % diffuse part of the Ashikhmin-Shirley model: every facet has its own
            % coefficient cd1, along its own normal,
            %   int (1 - (1 - n.v/2)^5) (n.v) v dw = 1573/2688 * pi * n
            const = orbitConst;
            d = au2km(1.0, const) * 10^3; % m
            sunB = [0.3, -0.2, 0.93]; sunB = sunB / norm(sunB);
            sat = testPortingFixes.threeFacets();
            NS = sat.normal * sunB';
            cd1 = 28/23 .* sat.Cd ./ pi .* (1 - sat.F0) .* (1 - (1 - NS / 2).^5);
            expected = sum(-const.S0 / const.c .* sat.area .* NS .* cd1 .* (1573/2688 * pi) .* sat.normal, 1);

            [~, cdAS] = srpAS(sat, sunB, d, const, 100);
            [~, cdUni] = srpASuni(sat, sunB, d, const, 100);
            testCase.verifyEqual(cdAS, expected, 'RelTol', 1e-12);
            testCase.verifyEqual(cdUni, expected, 'RelTol', 1e-12);

            % the same as the sum of the facets taken one by one
            total = zeros(1, 3);
            for i = 1:3
                [~, cdI] = srpAS(testPortingFixes.threeFacets(i), sunB, d, const, 100);
                total = total + cdI;
            end
            testCase.verifyEqual(cdAS, total, 'RelTol', 1e-12);
        end

        function srpASTiltedFacetsMatchUniformReference(testCase)
            % importance sampling (local frame of each facet) against the uniform
            % sampling in the body frame, for facets that are not normal to +z
            const = orbitConst;
            d = au2km(1.0, const) * 10^3; % m
            sunB = [0.3, -0.2, 0.93]; sunB = sunB / norm(sunB);
            rng(1, 'twister');
            satAS = srpAS(testPortingFixes.threeFacets(), sunB, d, const, 4e5);
            rng(1001, 'twister');
            satUni = srpASuni(testPortingFixes.threeFacets(), sunB, d, const, 1e6);
            err = vecnorm(satAS.force - satUni.force, 2, 2) ./ vecnorm(satUni.force, 2, 2);
            testCase.verifyLessThan(max(err), 0.02); % Monte-Carlo noise is about 0.1 %
        end

        function crossIsRowWiseForThreeRows(testCase)
            % cross(A, B) without dim works along the columns when A and B are 3 x 3,
            % i.e. for a model with exactly 3 facets or a history of exactly 3 steps
            sat = testPortingFixes.threeFacets();
            testCase.verifyEqual(calcNormalObj(sat), sat.normal, 'AbsTol', 1e-15);
            for i = 1:3
                testCase.verifyEqual(sat.uv(i,:), cross(sat.normal(i,:), sat.uu(i,:)), 'AbsTol', 1e-15);
            end

            const = orbitConst;
            d = au2km(1.0, const) * 10^3; % m
            sunB = [0.3, -0.2, 0.93]; sunB = sunB / norm(sunB);
            sat.Ca = 1 - sat.Cd - 0.2; sat.Cs = 0.2 * ones(3, 1); sat.kappa = 0; % for srpSimple
            sat.mCT = [0.1; 0.2; 0.3];                            % for srpCT, srpCTuni
            results = {srpAS(sat, sunB, d, const, 100), srpASuni(sat, sunB, d, const, 100), ...
                srpCT(sat, sunB, d, const, 'Beckmann', 100), ...
                srpCTuni(sat, sunB, d, const, 'Beckmann', 100), srpSimple(sat, sunB, d, const)};
            for k = 1:numel(results)
                s = results{k};
                for i = 1:3
                    testCase.verifyEqual(s.torque(i,:), cross(s.pos(i,:), s.force(i,:)), ...
                        'AbsTol', 1e-20);
                end
            end

            % 3 time steps: h = r x v of each row is perpendicular to r and v
            chief.oe = [7000, 1e-3, 0.9, 0.7, 0.4, 0.1; 7000, 1e-3, 0.9, 0.7, 0.4, 1.1; ...
                        7000, 1e-3, 0.9, 0.7, 0.4, 2.1];
            deputy.oe = chief.oe + [0.1, 1e-4, 1e-4, 2e-4, -3e-3, 2e-3];
            chief.n = sqrt(const.GE / 7000^3);
            chief.q = repmat([0 0 0 1], 3, 1); deputy.q = chief.q; deputy.w = zeros(3, 3);
            [chief4, ~, rel4] = calcRelPosVelAtti(testPortingFixes.addRow(chief), ...
                testPortingFixes.addRow(deputy), 1, const);
            [chief3, ~, rel3] = calcRelPosVelAtti(chief, deputy, 1, const);
            testCase.verifyEqual(chief3.qoi, chief4.qoi(1:3,:), 'AbsTol', 1e-12);
            testCase.verifyEqual(rel3.vNonlinRTN, rel4.vNonlinRTN(1:3,:), 'AbsTol', 1e-12);
        end
    end
end
