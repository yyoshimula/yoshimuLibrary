function dxdt = eomHCW(~, x_, n)

x = x_(:);
u = zeros(3,1);

A21 = [3*n^2, 0, 0
    0, 0, 0
    0, 0, -n^2];
A22 = [0, 2*n, 0
    -2*n, 0, 0
    0, 0, 0];

A = [zeros(3,3), eye(3)
    A21, A22];

B = [zeros(3,3)
    eye(3)];

dxdt = A * x + B * u;

end

