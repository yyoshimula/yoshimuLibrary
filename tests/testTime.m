classdef testTime < matlab.unittest.TestCase
    % Unit tests for time/ utilities.
    %
    % Reference values:
    %   - J2000.0 epoch: 2000-01-01 12:00:00 -> JD 2451545.0 (exact)
    %   - Meeus, Astronomical Algorithms, ch. 7 examples (JD table, p. 62)
    %   - Vallado, Fundamentals of Astrodynamics and Applications,
    %     Ex. 3-4 (1996-10-26 14:20 -> JD 2450383.09722222) and
    %     Ex. 3-7 leap seconds (mid-2004: deltaAT = 32 s)
    %   - tests/manual/verifyTimeCalc.m recorded SPICE output:
    %     2003-12-19 00:00:00 UTC -> JD(TT) 2452992.500742870
    %
    % Tolerances: JD values near 2.45e6 have ulp ~4.7e-10 day (~40 us), so
    % chained conversions are compared with AbsTol 1e-8 day (~1 ms); values
    % that are exact in binary use 1e-9 or exact comparison.

    properties (Constant)
        % official IERS TAI-UTC history since 1972
        % (https://hpiers.obspm.fr/iers/bul/bulc/TimeSteps.history):
        % [year month day] on which the new TAI-UTC takes effect (00:00 UTC,
        % i.e. the day after the leap-second insertion), and TAI-UTC [s]
        iersTaiUtc = [ ...
            1972  7  1  11
            1973  1  1  12
            1974  1  1  13
            1975  1  1  14
            1976  1  1  15
            1977  1  1  16
            1978  1  1  17
            1979  1  1  18
            1980  1  1  19
            1981  7  1  20
            1982  7  1  21
            1983  7  1  22
            1985  7  1  23
            1988  1  1  24
            1990  1  1  25
            1991  1  1  26
            1992  7  1  27
            1993  7  1  28
            1994  7  1  29
            1996  1  1  30
            1997  7  1  31
            1999  1  1  32
            2006  1  1  33
            2009  1  1  34
            2012  7  1  35
            2015  7  1  36
            2017  1  1  37];
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
    end

    methods (Test)

        function testJ2000Epoch(testCase)
            % J2000.0: 2000-01-01 12:00:00 -> JD 2451545.0
            testCase.verifyEqual(gc2jd(2000, 1, 1, 12, 0, 0), 2451545.0, ...
                'AbsTol', 1e-9);
        end

        function testKnownJulianDates(testCase)
            % Meeus ch.7 table values and Vallado Ex. 3-4
            testCase.verifyEqual(gc2jd(2000, 1, 1, 0, 0, 0), 2451544.5, ...
                'AbsTol', 1e-9);
            testCase.verifyEqual(gc2jd(1999, 1, 1, 0, 0, 0), 2451179.5, ...
                'AbsTol', 1e-9);
            testCase.verifyEqual(gc2jd(1987, 1, 27, 0, 0, 0), 2446822.5, ...
                'AbsTol', 1e-9);
            testCase.verifyEqual(gc2jd(1988, 6, 19, 12, 0, 0), 2447332.0, ...
                'AbsTol', 1e-9);
            % Meeus: Sputnik launch epoch, 1957 October 4.81 = JD 2436116.31
            testCase.verifyEqual(gc2jd(1957, 10, 4, 19, 26, 24), 2436116.31, ...
                'AbsTol', 1e-8);
            % Vallado Ex. 3-4: 1996-10-26 14:20:00 -> 2450383.09722222
            testCase.verifyEqual(gc2jd(1996, 10, 26, 14, 20, 0), ...
                2450383.09722222, 'AbsTol', 1e-8);
            % leap day (Jan/Feb branch of the Meeus formula)
            testCase.verifyEqual(gc2jd(2000, 2, 29, 0, 0, 0), 2451603.5, ...
                'AbsTol', 1e-9);
            testCase.verifyEqual(gc2jd(2000, 1, 31, 0, 0, 0), 2451574.5, ...
                'AbsTol', 1e-9);
        end

        function testGc2jdVectorized(testCase)
            % vector inputs give the same result as scalar calls
            y  = [2000; 1999; 1996];
            mo = [1; 1; 10];
            d  = [1; 1; 26];
            h  = [12; 0; 14];
            mi = [0; 0; 20];
            s  = [0; 0; 0];
            jdVec = gc2jd(y, mo, d, h, mi, s);
            jdScalar = zeros(3, 1);
            for k = 1:3
                jdScalar(k) = gc2jd(y(k), mo(k), d(k), h(k), mi(k), s(k));
            end
            testCase.verifyEqual(jdVec, jdScalar, 'AbsTol', 1e-9);
        end

        function testJd2gcJ2000(testCase)
            % JD 2451545.0 -> 2000-01-01 12:00:00
            [y, mo, d, h, mi, s] = jd2gc(2451545.0);
            testCase.verifyEqual(y, 2000, 'AbsTol', 0);
            testCase.verifyEqual(mo, 1, 'AbsTol', 0);
            testCase.verifyEqual(d, 1, 'AbsTol', 0);
            testCase.verifyEqual(h, 12, 'AbsTol', 0);
            testCase.verifyEqual(mi, 0, 'AbsTol', 0);
            testCase.verifyEqual(s, 0, 'AbsTol', 1e-6);
        end

        function testJdGcRoundTrip(testCase)
            % jd -> gregorian -> jd round trip to ~1 ms (a few ulp of JD)
            jdList = [2451545.0
                      2450383.0972222222
                      2436116.31
                      2447332.0
                      2446822.5
                      2452992.5];
            for k = 1:numel(jdList)
                [y, mo, d, h, mi, s] = jd2gc(jdList(k));
                jdBack = gc2jd(y, mo, d, h, mi, s);
                testCase.verifyEqual(jdBack, jdList(k), 'AbsTol', 1e-8, ...
                    sprintf('round trip failed for jd = %.10f', jdList(k)));
            end
        end

        function testMjdConversions(testCase)
            % MJD = JD - 2400000.5 (exact definition)
            testCase.verifyEqual(jd2mjd(2400000.5), 0, 'AbsTol', 1e-9);
            testCase.verifyEqual(jd2mjd(2451545.0), 51544.5, 'AbsTol', 1e-9);
            testCase.verifyEqual(mjd2jd(51544.5), 2451545.0, 'AbsTol', 1e-9);
            % round trip is exact in floating point
            jd = 2452992.500742870;
            testCase.verifyEqual(mjd2jd(jd2mjd(jd)), jd, 'AbsTol', 1e-9);
        end

        function testJulianCentury(testCase)
            % T = (jd - 2451545) / 36525: zero at J2000, one a Julian century later
            testCase.verifyEqual(jd2jdT(2451545.0), 0, 'AbsTol', 1e-12);
            testCase.verifyEqual(jd2jdT(2451545.0 + 36525), 1, 'AbsTol', 1e-12);
            testCase.verifyEqual(jd2jdT(2451545.0 - 36525), -1, 'AbsTol', 1e-12);
        end

        function testLeapSecondTableFormat(testCase)
            % leapS returns [JD from which the cumulative TAI-UTC applies
            % (00:00 UTC on the day AFTER the insertion), cumulative leap
            % seconds]; both columns must match the IERS history.
            leapJD = leapS();
            tbl = testTime.iersTaiUtc;
            nEntries = size(tbl, 1);
            z = zeros(nEntries, 1);
            expectedJD = gc2jd(tbl(:,1), tbl(:,2), tbl(:,3), z, z, z);
            testCase.verifySize(leapJD, [nEntries 2], ...
                'leapS must return one row per IERS leap second');
            testCase.verifyEqual(leapJD(:,1), expectedJD, 'AbsTol', 1e-9, ...
                'leap-second effective dates must match IERS');
            testCase.verifyEqual(leapJD(:,2) + 10, tbl(:,4), ...
                'cumulative leap seconds must match IERS TAI-UTC');
            % strictly increasing dates and cumulative seconds
            testCase.verifyGreaterThan(diff(leapJD(:, 1)), 0);
            testCase.verifyGreaterThan(diff(leapJD(:, 2)), 0);
        end

        function testDeltaAT(testCase)
            % deltaAT = TAI - UTC = 10 s + cumulative leap seconds
            leapJD = leapS();
            % before 1972: 10 s (by dAT convention)
            testCase.verifyEqual(dAT(gc2jd(1971, 6, 1, 0, 0, 0), leapJD), 10, ...
                'AbsTol', 0);
            % after the second leap second (1972-12-31): 12 s
            testCase.verifyEqual(dAT(gc2jd(1973, 1, 10, 0, 0, 0), leapJD), 12, ...
                'AbsTol', 0);
            % Vallado Ex. 3-7: May 2004 -> 32 s
            testCase.verifyEqual(dAT(gc2jd(2004, 5, 14, 0, 0, 0), leapJD), 32, ...
                'AbsTol', 0);
            % 2006-01-01 (after 2005-12-31 leap second): 33 s
            testCase.verifyEqual(dAT(gc2jd(2006, 1, 2, 0, 0, 0), leapJD), 33, ...
                'AbsTol', 0);
            % after 2017-01-01: 37 s
            testCase.verifyEqual(dAT(gc2jd(2017, 6, 1, 0, 0, 0), leapJD), 37, ...
                'AbsTol', 0);
            testCase.verifyEqual(dAT(gc2jd(2020, 1, 1, 0, 0, 0), leapJD), 37, ...
                'AbsTol', 0);
        end

        function testDeltaATAtAllIERSTransitions(testCase)
            % for every IERS entry: the previous TAI-UTC must hold through
            % the leap-second day itself (effective JD - 0.5 = 12:00 UTC on
            % the insertion day), and the new value from 00:00 UTC onward
            leapJD = leapS();
            tbl = testTime.iersTaiUtc;
            for k = 1:size(tbl, 1)
                effJD = gc2jd(tbl(k,1), tbl(k,2), tbl(k,3), 0, 0, 0);
                if k == 1
                    prevTaiUtc = 10;
                else
                    prevTaiUtc = tbl(k-1, 4);
                end
                tag = sprintf('entry %d (effective %04d-%02d-%02d)', ...
                    k, tbl(k,1), tbl(k,2), tbl(k,3));

                testCase.verifyEqual(dAT(effJD - 0.5, leapJD), ...
                    prevTaiUtc, ['noon of leap-second day, ' tag]);
                testCase.verifyEqual(dAT(effJD, leapJD), ...
                    tbl(k,4), ['start of effective day, ' tag]);
                testCase.verifyEqual(dAT(effJD + 0.5, leapJD), ...
                    tbl(k,4), ['noon of effective day, ' tag]);
            end
        end

        function testDeltaATBoundary2016to2017(testCase)
            % leap second inserted at 2016-12-31 23:59:60 UTC -> 37 s from
            % 2017-01-01 00:00 UTC, 36 s throughout 2016-12-31
            leapJD = leapS();
            testCase.verifyEqual(dAT(gc2jd(2016, 12, 30, 12, 0, 0), leapJD), 36, ...
                'day before the leap-second day');
            testCase.verifyEqual(dAT(gc2jd(2016, 12, 31, 12, 0, 0), leapJD), 36, ...
                'noon of the leap-second day');
            testCase.verifyEqual(dAT(gc2jd(2016, 12, 31, 23, 59, 59), leapJD), 36, ...
                'last second of the leap-second day');
            testCase.verifyEqual(dAT(gc2jd(2017, 1, 1, 0, 0, 0), leapJD), 37, ...
                'start of the day after');
            testCase.verifyEqual(dAT(gc2jd(2017, 1, 1, 12, 0, 0), leapJD), 37, ...
                'noon of the day after');
        end

        function testDeltaATBoundary2005to2006(testCase)
            % leap second inserted at 2005-12-31 23:59:60 UTC -> 33 s from
            % 2006-01-01 00:00 UTC, 32 s throughout 2005-12-31
            leapJD = leapS();
            testCase.verifyEqual(dAT(gc2jd(2005, 12, 31, 12, 0, 0), leapJD), 32, ...
                'noon of the leap-second day');
            testCase.verifyEqual(dAT(gc2jd(2005, 12, 31, 23, 59, 59), leapJD), 32, ...
                'last second of the leap-second day');
            testCase.verifyEqual(dAT(gc2jd(2006, 1, 1, 0, 0, 0), leapJD), 33, ...
                'start of the day after');
        end

        function testDeltaATVectorInput(testCase)
            % dAT must accept a jd column vector (documented jd (:,1)) and
            % return elementwise results matching the scalar calls
            leapJD = leapS();
            jd = [gc2jd(1971,  6, 30, 12, 0, 0)
                  gc2jd(2005, 12, 31, 12, 0, 0)
                  gc2jd(2006,  1,  1,  0, 0, 0)
                  gc2jd(2016, 12, 31, 12, 0, 0)
                  gc2jd(2017,  1,  1,  0, 0, 0)
                  gc2jd(2020,  1,  1, 12, 0, 0)];
            expected = [10; 32; 33; 36; 37; 37];

            actual = dAT(jd, leapJD);
            testCase.verifySize(actual, size(jd), ...
                'deltaAT must have the same size as jd');
            testCase.verifyEqual(actual, expected);

            for k = 1:numel(jd)
                testCase.verifyEqual(dAT(jd(k), leapJD), expected(k), ...
                    'vector result must match the scalar call');
            end
        end

        function testUtc2ttOffsetRelation(testCase)
            % TT = UTC + deltaAT + 32.184 s.
            % AbsTol 2e-4 s: the JD sum rounds at ulp(2.45e6 days) ~ 40 us per
            % operation, so the recovered offset is good to ~1e-4 s.
            leapJD = leapS();
            dates = [2004, 5, 14, 10, 43, 0
                     2017, 6, 1, 0, 0, 0
                     1995, 3, 30, 6, 15, 30];
            for k = 1:size(dates, 1)
                jdUTC = gc2jd(dates(k, 1), dates(k, 2), dates(k, 3), ...
                    dates(k, 4), dates(k, 5), dates(k, 6));
                deltaAT = dAT(jdUTC, leapJD);
                jdTT = utc2tt(jdUTC, deltaAT);
                offsetSec = (jdTT - jdUTC) * 86400;
                testCase.verifyEqual(offsetSec, deltaAT + 32.184, ...
                    'AbsTol', 2e-4, ...
                    sprintf('TT-UTC offset wrong for row %d', k));
            end
        end

        function testUtc2ttSpiceReference(testCase)
            % 2003-12-19 00:00:00 UTC -> JD(TT) 2452992.500742870
            % (SPICE JDTDT value recorded in tests/manual/verifyTimeCalc.m;
            % deltaAT = 32 s at that epoch)
            leapJD = leapS();
            jdUTC = gc2jd(2003, 12, 19, 0, 0, 0);
            deltaAT = dAT(jdUTC, leapJD);
            testCase.verifyEqual(deltaAT, 32, 'AbsTol', 0);
            jdTT = utc2tt(jdUTC, deltaAT);
            testCase.verifyEqual(jdTT, 2452992.500742870, 'AbsTol', 1e-8);
        end

        function testUt2ttTableValue(testCase)
            % at J2000.0 the interpolation lands exactly on the year-2000 table
            % node (n = 0), so ut2tt returns the tabulated 63.8285 s exactly
            testCase.verifyEqual(ut2tt(2451545.0), 63.8285, 'AbsTol', 1e-9);
        end

        function testDoy2gc(testCase)
            % 2004 day-of-year 32.5 -> Feb 1, 12:00:00 (Jan has 31 days)
            [mo, d, h, mi, s] = doy2gc(2004, 32.5);
            testCase.verifyEqual(mo, 2, 'AbsTol', 0);
            testCase.verifyEqual(d, 1, 'AbsTol', 0);
            testCase.verifyEqual(h, 12, 'AbsTol', 0);
            testCase.verifyEqual(mi, 0, 'AbsTol', 0);
            testCase.verifyEqual(s, 0, 'AbsTol', 1e-4);
            % doy 1.0 is Jan 1 00:00
            [mo1, d1, h1, mi1, s1] = doy2gc(2004, 1);
            testCase.verifyEqual([mo1, d1, h1, mi1], [1, 1, 0, 0], 'AbsTol', 0);
            testCase.verifyEqual(s1, 0, 'AbsTol', 1e-4);
        end

    end
end
