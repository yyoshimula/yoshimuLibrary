classdef testSRP < matlab.unittest.TestCase
    % testSRP  Regression and smoke tests for the SRP-related library code.
    %
    % Regression baselines for ctM were captured from the pre-consolidation
    % implementation (srp/ctM.m before it was merged into hifiSRP/ctM.m).
    % They are ground truth: if a ctM test fails, the consolidation broke
    % something. Do NOT loosen tolerances and do NOT edit ctM.m -- report it.
    %
    % Run headless with:
    %   matlab -batch "r = runtests('tests/testSRP.m'); disp(table(r)); assertSuccess(r)"

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
        function [sat, v, sunB] = singleFacetCase()
            % Common single-facet configuration for the ctM regression tests.
            sat.normal = [0.0, sin(0.3), cos(0.3)];
            sat.area = 1;
            sat.F0 = 0.5;
            sat.mCT = 0.05;

            v = [0.3 -0.2 0.933];
            v = v / norm(v);
            sunB = [0.1 0.4 0.911];
            sunB = sunB / norm(sunB);
        end

        function sat = interpTestSatellite()
            % Single flat facet used by tests/manual/verifyInterp.m.
            nFacet = 1;
            sat.pos = [0 0 0];
            sat.normal = [0 0 1];
            sat.area = 1;
            sat.qlb = [0 0 0 1];
            sat.F0 = 0.5 .* ones(nFacet, 1);
            sat.kappa = 0 .* ones(nFacet, 1);
            sat.Cd = 0.5 .* ones(nFacet, 1);
            sat.Cs = sat.F0;
            sat.Ca = 0.0 .* ones(nFacet, 1);
            sat.mCT = 0.05 .* ones(nFacet, 1);
        end
    end

    methods (Test)
        function testCtMSingleFacetScalarF0(testCase)
            % (a) ctM regression, single facet, scalar F0.
            [sat, v, sunB] = testSRP.singleFacetCase();

            M = ctM(sat, v, sunB);

            expected = 0.12495897451119992; % pre-consolidation baseline
            testCase.verifyEqual(M, expected, 'RelTol', 1e-12, ...
                'ctM single-facet regression value changed after consolidation.');
        end

        function testCtMSingleFacetVectorF0(testCase)
            % (b) Same as (a) but F0 given as a 1x1 array.
            [sat, v, sunB] = testSRP.singleFacetCase();
            sat.F0 = 0.5 * ones(1, 1);

            M = ctM(sat, v, sunB);

            expected = 0.12495897451119992; % pre-consolidation baseline
            testCase.verifyEqual(M, expected, 'RelTol', 1e-12, ...
                'ctM must give the same result for scalar and 1x1 F0.');
        end

        function testCtMMultiFacetScalarF0(testCase)
            % (c) ctM regression, multiple facets with a scalar F0.
            [~, v, sunB] = testSRP.singleFacetCase();

            satM.normal = [0 0 1; 0 1 0; 1 0 0; 0 0 -1];
            satM.area = ones(4, 1);
            satM.F0 = 0.5;
            satM.mCT = 0.05;

            M = ctM(satM, v, sunB);

            testCase.verifySize(M, [1 1], ...
                'ctM must sum over facets and return a 1x1 scalar.');
            expected = 0.24991676671001906; % pre-consolidation baseline
            testCase.verifyEqual(M, expected, 'RelTol', 1e-12, ...
                'ctM multi-facet regression value changed after consolidation.');
        end

        function testCtMAcceptsRowOrColumnVectors(testCase)
            % (d) ctM must accept both row and column v / sunB.
            [sat, v, sunB] = testSRP.singleFacetCase();

            Mrow = ctM(sat, v, sunB);
            Mcol = ctM(sat, v', sunB');

            testCase.verifyEqual(Mcol, Mrow, 'RelTol', 1e-14, ...
                'ctM must be insensitive to row/column orientation of v and sunB.');
        end

        function testSrpCTinterpMatchesSrpCT(testCase)
            % (e) Smoke test: srpCTinterp (correction-table interpolation)
            % must loosely reproduce the Monte-Carlo srpCT specular part.
            % Mirrors tests/manual/verifyInterp.m (which compares the two and
            % plots absolute/relative error); here we only assert a loose
            % agreement: 30 % relative per nonzero component, absolute 1e-9
            % where the reference component is essentially zero (the x
            % component is ~1e-10 N for phiI = -90 deg). This checks that the
            % interpolation table matches its consumer, not precision.
            s = load('correctionPara.mat'); % lives at srp/correctionPara.mat, found via path
            correctionPara = s.correctionPara;

            sat = testSRP.interpTestSatellite();

            const = orbitConst;
            d = au2km(1.0, const) * 10^3; % m, distance from sat to sun

            phiI = deg2rad(-90);
            absTol = 1e-9; % N, for near-zero reference components
            relTol = 0.30;

            for thetaIdeg = [10 30 60]
                thetaI = deg2rad(thetaIdeg);
                sunB = [sin(thetaI)*cos(phiI), sin(thetaI)*sin(phiI), cos(thetaI)];

                % srpCT integrates by Monte Carlo; fix the seed so the
                % reference value is deterministic (no user input, no figures).
                rng(12345, 'twister');
                [~, ~, srpCsCT] = srpCT(sat, sunB, d, const);
                [~, ~, srpCsApprox] = srpCTinterp(sat, sunB, d, const, correctionPara);

                testCase.verifySize(srpCsCT, [1 3], ...
                    sprintf('srpCT specular output must be 1x3 (thetaI = %d deg).', thetaIdeg));
                testCase.verifySize(srpCsApprox, [1 3], ...
                    sprintf('srpCTinterp specular output must be 1x3 (thetaI = %d deg).', thetaIdeg));
                testCase.verifyTrue(all(isfinite(srpCsCT)), ...
                    sprintf('srpCT returned non-finite values (thetaI = %d deg).', thetaIdeg));
                testCase.verifyTrue(all(isfinite(srpCsApprox)), ...
                    sprintf('srpCTinterp returned non-finite values (thetaI = %d deg).', thetaIdeg));

                for k = 1:3
                    diagInfo = sprintf(['srpCTinterp vs srpCT specular component %d ' ...
                        'at thetaI = %d deg (ref = %.6e, interp = %.6e).'], ...
                        k, thetaIdeg, srpCsCT(k), srpCsApprox(k));
                    if abs(srpCsCT(k)) > absTol
                        relErr = abs(srpCsApprox(k) - srpCsCT(k)) / abs(srpCsCT(k));
                        testCase.verifyLessThan(relErr, relTol, diagInfo);
                    else
                        testCase.verifyLessThan(abs(srpCsApprox(k) - srpCsCT(k)), absTol, diagInfo);
                    end
                end
            end
        end

        function testLibraryPathNotShadowed(testCase)
            % (g) Guard for addLibraryToPath: stale library copies under
            % hidden folders (e.g. .claude/worktrees session worktrees)
            % must not shadow the real library files, otherwise the whole
            % suite silently tests outdated code.
            root = fileparts(fileparts(mfilename('fullpath')));
            testCase.verifyEqual(which('srpCT'), fullfile(root, 'srp', 'srpCT.m'), ...
                'srpCT resolves to a shadowed copy; check addLibraryToPath.');
        end

        function testSrpCTGaussMatchesUniformReference(testCase)
            % (h) Cross-validation of the srpCT Gauss NDF branch against
            % the uniform-sampling reference srpCTuni (same Cook-Torrance
            % model, Gaussian NDF), following the methodology of
            % tests/manual/verifyImpSampling.m. Guards the half-vector
            % normalization order and the branch-specific Monte-Carlo
            % weight (unlike the Beckmann branch, D and the sin(thetaR)
            % Jacobian do not cancel for uniform sampling): before the
            % 2026-07 fix the two disagreed by >100 % with opposite signs
            % in x and y. Both estimates are seeded, so the comparison is
            % deterministic; 10 % covers the Monte-Carlo noise (observed
            % disagreement is below 3 %).
            sat = testSRP.interpTestSatellite();
            sat.mCT = 0.3; % broad specular lobe so both MC estimates converge

            const = orbitConst;
            d = au2km(1.0, const) * 10^3; % m, distance from sat to sun

            thetaI = deg2rad(35);
            phiI = deg2rad(20);
            sunB = [sin(thetaI)*cos(phiI), sin(thetaI)*sin(phiI), cos(thetaI)];

            rng(1, 'twister');
            [~, ~, csGauss] = srpCT(sat, sunB, d, const, 'Gauss', 2e5);
            rng(1001, 'twister');
            [~, ~, csUni] = srpCTuni(sat, sunB, d, const, 'Gauss', 4e5);

            testCase.verifySize(csGauss, [1 3], ...
                'srpCT Gauss specular output must be 1x3.');
            testCase.verifyTrue(all(isfinite(csGauss)), ...
                'srpCT Gauss branch returned non-finite values.');
            testCase.verifyEqual(csGauss, csUni, 'RelTol', 0.10, ...
                'srpCT Gauss branch disagrees with the srpCTuni uniform-sampling reference.');
        end

        function testSrpCTGaussRegression(testCase)
            % (i) Deterministic regression pin for the srpCT Gauss branch
            % (seeded Monte Carlo). Baseline captured from the fixed
            % implementation (2026-07, after the half-vector normalization
            % and MC-weight corrections). If an intentional algorithm
            % change moves these values, re-baseline ONLY after
            % re-verifying against srpCTuni (test h).
            sat = testSRP.interpTestSatellite();

            const = orbitConst;
            d = au2km(1.0, const) * 10^3; % m, distance from sat to sun

            thetaI = deg2rad(30);
            phiI = deg2rad(-90);
            sunB = [sin(thetaI)*cos(phiI), sin(thetaI)*sin(phiI), cos(thetaI)];

            rng(12345, 'twister');
            [~, srpCdOut, srpCsOut] = srpCT(sat, sunB, d, const, 'Gauss', 1e4);

            expectedCs = [-1.0358329333471533e-10, -6.9932828453404173e-09, -1.2237565252850401e-08];
            testCase.verifyEqual(srpCsOut, expectedCs, 'RelTol', 1e-10, ...
                'Gauss-branch specular regression value changed.');

            % diffuse part is analytic and branch-independent
            expectedCd = [0, 0, -1.3066778272492802e-06];
            testCase.verifyEqual(srpCdOut, expectedCd, 'RelTol', 1e-10, ...
                'AbsTol', 1e-20, 'Diffuse SRP regression value changed.');
        end

        function testMagDimsWithDistance(testCase)
            % (f-1) lightcurves/mag.m: apparent magnitude must increase
            % (object gets dimmer) as distance d grows with fObs fixed.
            fObs = 1e-10;
            dSpan = [1e5, 1e6, 1e7, 1e8]; % m

            appMag = mag(fObs, dSpan);

            testCase.verifyTrue(all(diff(appMag) > 0), ...
                'Apparent magnitude must increase (dim) with distance for fixed fObs.');
        end

        function testMagHandComputedValue(testCase)
            % (f-2) Hand-computed value from mag.m's formula:
            % appMag = -26.7 - 2.5*log10(fObs / d^2)
            % with fObs = 1e-12, d = 1e3 m:
            %   fObs/d^2 = 1e-18, log10 = -18, appMag = -26.7 + 45 = 18.3.
            appMag = mag(1e-12, 1e3);

            testCase.verifyEqual(appMag, 18.3, 'AbsTol', 1e-10, ...
                'mag.m formula regression value changed.');
        end
    end
end
