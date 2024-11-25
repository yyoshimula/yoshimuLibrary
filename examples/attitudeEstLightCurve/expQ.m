function q = expQ(scalar, x)
arguments
    scalar (1,1) {mustBeMember(scalar, [0, 4])}
    x (:,3) {mustBeNumeric}
end

xNorm = vecnorm(x, 2, 2);
q = (scalar == 0) .* [cos(xNorm), normRow(x) .* sin(xNorm)] ...
    + (scalar == 4) .* [normRow(x) .* sin(xNorm), cos(xNorm)];

end

