function dxdt = eulerEom(t, x, MOI)
% -------------------------------------------------------------------------------------
%   Equations of motion for rigid body
%   12ŒŽ 2, 2013 yasuhiro yoshimura
%    Input:
%   Output:
%   related function files: used in eulerMain.m
%-------------------------------------------------------------------------------------


% state variables
q1 = x(1);
q2 = x(2);
q3 = x(3);
q4 = x(4);
q = [q1;q2;q3;q4];

wx = x(5);
wy = x(6);
wz = x(7);

wVec = [wx; wy; wz];

%------Equations of Motion---------------------

%--Quaternion-----------------

% control inputs
% u = -[q1;q2;q3] - wVec;

% u = (- 4.*wVec + 2.*[q1; q2; q3]);
% u = u./ norm(u) .* 1e-4;
u = zeros(3,1);

%------->state variables [q0, q1, q2, q3, wx, wy, wz]
dxdt = [ 0.5.*qMult(4, 1, [wVec;0]', q')'
    -MOI^(-1) * cross(wVec, (MOI * wVec)) + MOI^(-1) * u];