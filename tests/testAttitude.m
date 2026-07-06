classdef testAttitude < matlab.unittest.TestCase
    % Unit tests for attitude/ and math/ utilities.
    %
    % Conventions under test (see attitude/README.md and function headers):
    %   scalar == 0 : q = [q0, q1, q2, q3] = [cos(th/2), e*sin(th/2)] (scalar first)
    %   scalar == 4 : q = [q1, q2, q3, q4] = [e*sin(th/2), cos(th/2)] (scalar last)
    %   q2dcm returns the frame-rotation (attitude) matrix, i.e. for a rotation
    %   by th about x: A = [1 0 0; 0 c s; 0 -s c] = dcm1axis(1, th).
    %   qMult def == 1 (otimes): A(q ox p) = A(q) * A(p)
    %   qMult def == 0 (odot)  : A(q od p) = A(p) * A(q)
    %
    % Tolerances: DCM entries are short products of O(1) doubles, so a few eps
    % (~2e-16) of error per entry is expected; 1e-15..1e-14 covers 5-50 eps.
    % Chained conversions (round trips through trig) get 1e-12.

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

        function testQ2dcmOrthonormalRandom(testCase)
            % DCM from a random unit quaternion is orthonormal with det = +1
            rng(20260706, 'twister');
            for k = 1:20
                q = randn(1, 4);
                q = q ./ norm(q);
                for scalar = [0 4]
                    R = q2dcm(scalar, q);
                    testCase.verifyEqual(R' * R, eye(3), 'AbsTol', 1e-14, ...
                        sprintf('R''*R ~= I for scalar=%d, k=%d', scalar, k));
                    testCase.verifyEqual(det(R), 1.0, 'AbsTol', 1e-14, ...
                        sprintf('det(R) ~= +1 for scalar=%d, k=%d', scalar, k));
                end
            end
        end

        function testQ2dcmIdentityQuaternion(testCase)
            % identity quaternion -> identity matrix
            testCase.verifyEqual(q2dcm(0, [1 0 0 0]), eye(3), 'AbsTol', 1e-15);
            testCase.verifyEqual(q2dcm(4, [0 0 0 1]), eye(3), 'AbsTol', 1e-15);
        end

        function testQ2dcm90DegRotations(testCase)
            % 90 deg frame rotations about x, y, z (half angle = 45 deg)
            c = cos(pi/4);
            s = sin(pi/4);

            Rx90 = [1 0 0; 0 0 1; 0 -1 0];
            Ry90 = [0 0 -1; 0 1 0; 1 0 0];
            Rz90 = [0 1 0; -1 0 0; 0 0 1];

            % scalar-first convention
            testCase.verifyEqual(q2dcm(0, [c s 0 0]), Rx90, 'AbsTol', 1e-15);
            testCase.verifyEqual(q2dcm(0, [c 0 s 0]), Ry90, 'AbsTol', 1e-15);
            testCase.verifyEqual(q2dcm(0, [c 0 0 s]), Rz90, 'AbsTol', 1e-15);

            % scalar-last convention
            testCase.verifyEqual(q2dcm(4, [s 0 0 c]), Rx90, 'AbsTol', 1e-15);
            testCase.verifyEqual(q2dcm(4, [0 s 0 c]), Ry90, 'AbsTol', 1e-15);
            testCase.verifyEqual(q2dcm(4, [0 0 s c]), Rz90, 'AbsTol', 1e-15);
        end

        function testQ2dcmMatchesDcm1axis(testCase)
            % q1axis -> q2dcm equals dcm1axis / dcm1axisX/Y/Z for single-axis rotations
            rng(20260706, 'twister');
            axes3 = eye(3);
            singleAxisDcm = {@dcm1axisX, @dcm1axisY, @dcm1axisZ};
            for k = 1:5
                th = (2 * rand - 1) * pi;
                for ax = 1:3
                    Rref = dcm1axis(ax, th);
                    testCase.verifyEqual(singleAxisDcm{ax}(th), Rref, 'AbsTol', 1e-15);
                    for scalar = [0 4]
                        q = q1axis(scalar, axes3(ax, :), th);
                        testCase.verifyEqual(q2dcm(scalar, q), Rref, 'AbsTol', 1e-14, ...
                            sprintf('axis=%d, scalar=%d, th=%g', ax, scalar, th));
                    end
                end
            end
        end

        function testScalarFirstLastConsistency(testCase)
            % same physical rotation expressed in both conventions -> same DCM
            rng(20260706, 'twister');
            for k = 1:20
                q0first = randn(1, 4);
                q0first = q0first ./ norm(q0first);
                q4last = [q0first(2:4), q0first(1)];
                testCase.verifyEqual(q2dcm(0, q0first), q2dcm(4, q4last), ...
                    'AbsTol', 1e-15);
            end
        end

        function testDcm2qRoundTrip(testCase)
            % q -> DCM -> q recovers the quaternion up to the q <-> -q sign
            rng(20260706, 'twister');
            for k = 1:20
                q = randn(1, 4);
                q = q ./ norm(q);
                for scalar = [0 4]
                    if scalar == 0
                        qs = q;
                    else
                        qs = [q(2:4), q(1)];
                    end
                    qBack = dcm2q(scalar, q2dcm(scalar, qs));
                    testCase.verifySize(qBack, [1 4]);
                    errUpToSign = min(norm(qBack - qs), norm(qBack + qs));
                    testCase.verifyEqual(errUpToSign, 0, 'AbsTol', 1e-12, ...
                        sprintf('dcm2q round trip failed, scalar=%d, k=%d', scalar, k));
                end
            end
        end

        function testSkew(testCase)
            % skew(a)*b == cross(a,b); skew(a) is antisymmetric with zero diagonal
            rng(20260706, 'twister');
            for k = 1:10
                a = randn(3, 1);
                b = randn(3, 1);
                S = skew(a);
                testCase.verifyEqual(S * b, cross(a, b), 'AbsTol', 1e-15);
                testCase.verifyEqual(S', -S, 'AbsTol', 0);
                testCase.verifyEqual(diag(S), zeros(3, 1), 'AbsTol', 0);
            end
        end

        function testQMultMatMatchesQMult(testCase)
            % qMultMat(scalar, def, q) * p' == qMult(scalar, def, q, p)'
            rng(20260706, 'twister');
            for k = 1:10
                q = randn(1, 4);
                q = q ./ norm(q);
                p = randn(1, 4);
                p = p ./ norm(p);
                for scalar = [0 4]
                    if scalar == 0
                        qs = q;  ps = p;
                    else
                        qs = [q(2:4), q(1)];
                        ps = [p(2:4), p(1)];
                    end
                    for def = [0 1]
                        lhs = qMultMat(scalar, def, qs) * ps';
                        rhs = qMult(scalar, def, qs, ps)';
                        testCase.verifyEqual(lhs, rhs, 'AbsTol', 1e-15, ...
                            sprintf('scalar=%d, def=%d, k=%d', scalar, def, k));
                    end
                end
            end
        end

        function testQMultComposition(testCase)
            % def==1 (otimes): A(q ox p) = A(q)A(p)
            % def==0 (odot)  : A(q od p) = A(p)A(q)
            rng(20260706, 'twister');
            for k = 1:10
                q = randn(1, 4);
                q = q ./ norm(q);
                p = randn(1, 4);
                p = p ./ norm(p);
                for scalar = [0 4]
                    if scalar == 0
                        qs = q;  ps = p;
                    else
                        qs = [q(2:4), q(1)];
                        ps = [p(2:4), p(1)];
                    end
                    Aq = q2dcm(scalar, qs);
                    Ap = q2dcm(scalar, ps);
                    testCase.verifyEqual(q2dcm(scalar, qMult(scalar, 1, qs, ps)), ...
                        Aq * Ap, 'AbsTol', 1e-14, ...
                        sprintf('otimes composition, scalar=%d, k=%d', scalar, k));
                    testCase.verifyEqual(q2dcm(scalar, qMult(scalar, 0, qs, ps)), ...
                        Ap * Aq, 'AbsTol', 1e-14, ...
                        sprintf('odot composition, scalar=%d, k=%d', scalar, k));
                end
            end
        end

        function testQConjQInvUnitQuaternion(testCase)
            % for unit quaternions qInv == qConj and q * qInv == identity
            rng(20260706, 'twister');
            idQ = {[1 0 0 0], [0 0 0 1]};   % identity, scalar-first / scalar-last
            scalars = [0 4];
            for k = 1:10
                q = randn(1, 4);
                q = q ./ norm(q);
                for is = 1:2
                    scalar = scalars(is);
                    if scalar == 0
                        qs = q;
                    else
                        qs = [q(2:4), q(1)];
                    end
                    qi = qInv(scalar, qs);
                    testCase.verifyEqual(qi, qConj(scalar, qs), 'AbsTol', 1e-15);
                    for def = [0 1]
                        testCase.verifyEqual(qMult(scalar, def, qs, qi), ...
                            idQ{is}, 'AbsTol', 1e-15, ...
                            sprintf('q*qInv ~= identity, scalar=%d, def=%d', scalar, def));
                    end
                end
            end
        end

        function testQRotationMatchesDcm(testCase)
            % qRotation expresses a vector in the rotated frame: rb == A(q) * r
            rng(20260706, 'twister');
            for k = 1:10
                q = randn(1, 4);
                q = q ./ norm(q);
                r = randn(1, 3);
                for scalar = [0 4]
                    if scalar == 0
                        qs = q;
                    else
                        qs = [q(2:4), q(1)];
                    end
                    rb = qRotation(scalar, r, qs);
                    testCase.verifyEqual(rb', q2dcm(scalar, qs) * r', ...
                        'AbsTol', 1e-14, sprintf('scalar=%d, k=%d', scalar, k));
                end
            end
        end

        function testZyxEulerRoundTrip(testCase)
            % zyx2q -> q2zyx recovers angles inside the principal ranges;
            % zyx2dcm agrees with q2dcm(zyx2q) and with Rx(psi)*Ry(th)*Rz(phi)
            rng(20260706, 'twister');
            for k = 1:15
                phi   = (2 * rand - 1) * (pi - 0.1);      % 1st rot, about z
                theta = (2 * rand - 1) * (pi/2 - 0.1);    % 2nd rot, about y
                psi   = (2 * rand - 1) * (pi - 0.1);      % 3rd rot, about x
                Rref = dcm1axisX(psi) * dcm1axisY(theta) * dcm1axisZ(phi);
                testCase.verifyEqual(zyx2dcm(phi, theta, psi), Rref, ...
                    'AbsTol', 1e-14);
                for scalar = [0 4]
                    q = zyx2q(scalar, phi, theta, psi);
                    testCase.verifyEqual(q2dcm(scalar, q), Rref, 'AbsTol', 1e-14, ...
                        sprintf('zyx2q vs dcm, scalar=%d, k=%d', scalar, k));
                    testCase.verifyEqual(q2zyx(scalar, q), [phi, theta, psi], ...
                        'AbsTol', 1e-12, ...
                        sprintf('euler round trip, scalar=%d, k=%d', scalar, k));
                end
            end
        end

        function testRotVecRoundTrip(testCase)
            % rotVec2q -> q2rotVec recovers the rotation vector for angles in (0, pi)
            rng(20260706, 'twister');
            for k = 1:10
                e = randn(1, 3);
                e = e ./ norm(e);
                th = 0.2 + 2.7 * rand;   % in (0.2, 2.9) subset of (0, pi)
                rv = th .* e;
                for scalar = [0 4]
                    rvBack = q2rotVec(scalar, rotVec2q(scalar, rv));
                    testCase.verifyEqual(rvBack, rv, 'AbsTol', 1e-12, ...
                        sprintf('scalar=%d, k=%d', scalar, k));
                end
            end
        end

    end
end
