classdef testFrames < matlab.unittest.TestCase
    % Unit tests for the ITRF -> GCRF transformation (CIO-based
    % IAU-2006/2000 reduction) and related time/angle utilities.
    %
    % Reference case: Vallado, "Fundamentals of Astrodynamics and
    % Applications", 4th ed., Example 3-14 (see also
    % tests/manual/verifyITRF2GCRF.m):
    %   epoch  2004-04-06 07:51:28.386009 UTC
    %   rITRF = [-1033.4793830; 7901.2952754; 6380.3565958] km
    %   rGCRF = [ 5102.508959;  6123.011403;  6378.136925] km
    %
    % EOP data: orbit/EOP_20_C04_one_file_1962-now.txt (on the path).
    % The EOP file is parsed once in TestClassSetup and cached.

    properties (Constant)
        rITRF = [-1033.4793830; 7901.2952754; 6380.3565958]; % km
        rGCRFvallado = [5102.508959; 6123.011403; 6378.136925]; % km
    end

    properties
        EOP % cached output of readEOP (dataAll, iau06, leapJD)
        jdUTC % Julian day of the Vallado Example 3-14 epoch (UTC)
        dcm % DCM from ITRF to GCRF at the test epoch
    end

    methods (TestClassSetup)
        function addLibraryToPath(testCase) %#ok<MANU>
            root = fileparts(fileparts(mfilename('fullpath')));
            % genpath includes hidden folders such as .claude/worktrees,
            % whose stale library copies would shadow the real files
            p = strsplit(genpath(root), pathsep);
            p = p(~cellfun(@isempty, p) & ~contains(p, [filesep '.']));
            addpath(strjoin(p, pathsep));
        end

        function loadEopAndDcm(testCase)
            % Path setup is repeated here (idempotent) so this method does
            % not depend on the execution order of TestClassSetup methods.
            root = fileparts(fileparts(mfilename('fullpath')));
            % genpath includes hidden folders such as .claude/worktrees,
            % whose stale library copies would shadow the real files
            p = strsplit(genpath(root), pathsep);
            p = p(~cellfun(@isempty, p) & ~contains(p, [filesep '.']));
            addpath(strjoin(p, pathsep));

            testCase.EOP = readEOP("EOP_20_C04_one_file_1962-now.txt");
            testCase.jdUTC = gc2jd(2004, 4, 6, 7, 51, 28.386009);
            testCase.dcm = itrf2gcrf(testCase.jdUTC, testCase.EOP);
        end
    end

    methods (Test)
        %% ITRF -> GCRF, Vallado Example 3-14 -----------------------------
        function testItrf2GcrfValladoExample(testCase)
            % The library currently reproduces Vallado's rGCRF to ~6e-6 km;
            % 1e-3 km leaves margin for future EOP data-file updates.
            rGCRF = testCase.dcm * testCase.rITRF;
            testCase.verifyEqual(rGCRF, testCase.rGCRFvallado, ...
                'AbsTol', 1e-3, ...
                'ITRF->GCRF disagrees with Vallado Example 3-14');
        end

        function testDcmIsProperOrthonormal(testCase)
            % A rotation matrix satisfies R'*R = I and det(R) = +1.
            R = testCase.dcm;
            testCase.verifyEqual(R'*R, eye(3), 'AbsTol', 1e-12);
            testCase.verifyEqual(det(R), 1, 'AbsTol', 1e-12);
        end

        function testNormPreservation(testCase)
            % A frame rotation must preserve the vector norm.
            rGCRF = testCase.dcm * testCase.rITRF;
            testCase.verifyEqual(norm(rGCRF), norm(testCase.rITRF), ...
                'AbsTol', 1e-9);
        end

        %% time utilities used by the reduction ---------------------------
        function testGc2jdValladoEpoch(testCase)
            % JD(UTC) of the Example 3-14 epoch:
            % 2453101.5 + (7 + 51/60 + 28.386009/3600)/24
            jdExpected = 2453101.5 + (7 + 51/60 + 28.386009/3600)/24;
            testCase.verifyEqual(testCase.jdUTC, jdExpected, 'AbsTol', 1e-9);
        end

        function testDeltaATAtEpoch(testCase)
            % Vallado Example 3-14: deltaAT = TAI - UTC = 32 s in April 2004.
            deltaAT = dAT(testCase.jdUTC, testCase.EOP.leapJD);
            testCase.verifyEqual(deltaAT, 32, 'AbsTol', 0);
        end

        %% era / gmst sanity checks ---------------------------------------
        function testEraAtJ2000(testCase)
            % IAU-2000 Earth rotation angle at the J2000.0 epoch (UT1):
            % ERA(2451545.0) = 2*pi*0.7790572732640 rad (IERS Conventions,
            % also Vallado 4th ed., Eq. 3-55 with Tu = 0).
            theta = era(2451545.0);
            thetaExpected = mod(2*pi*0.7790572732640, 2*pi);
            testCase.verifyEqual(theta, thetaExpected, 'AbsTol', 1e-12);
        end

        function testGmstMeeusExample(testCase)
            % Meeus, "Astronomical Algorithms" 2nd ed., Example 12.b:
            % 1987 April 10, 19:21:00 UT -> GMST = 8h 34m 57.0896s.
            % The reference is quoted to 0.0001 s of time (= 7.3e-9 rad),
            % so 1e-8 rad is the tightest justified tolerance.
            jd = gc2jd(1987, 4, 10, 19, 21, 0);
            gmstExpected = (8 + 34/60 + 57.0896/3600)/24 * 2*pi; % rad
            testCase.verifyEqual(gmst(jd), gmstExpected, 'AbsTol', 1e-8);
        end
    end
end
