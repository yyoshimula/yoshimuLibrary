function dxdt = eulerEom(~, x, MOI)
% ----------------------------------------------------------------------
%   Equations of motion for rigid body
%    20240909  y.yoshimura
%    Inputs:
%   Outputs:
%   related function files:
%   note:
%   cf:
%   revisions;
%
%   (c) 2024 yasuhiro yoshimura
%----------------------------------------------------------------------

%% state variables
q1 = x(1);
q2 = x(2);
q3 = x(3);
q4 = x(4);

q = [q1, q2, q3, q4];

wx = x(5);
wy = x(6);
wz = x(7);


w = [wx, wy, wz];

%%------Equations of Motion---------------------

% control inputs
u = [0.0
    0.0
    0.0];

dxdt = [ qKine(4, q, w)'
    -MOI^(-1) * cross(w', (MOI * w')) + MOI^(-1) * u];

end