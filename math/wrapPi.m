function lambda = wrapPi(lambda)
% WRAPPI Wrap angle in radians to the interval [-pi, pi)
%
%   lambdaWrapped = wrapPi(LAMBDA) wraps angles in LAMBDA, in radians,
%   to the interval [-pi, pi).
%
%   See also mod.

lambda = mod(lambda + pi, 2*pi) - pi;
end
