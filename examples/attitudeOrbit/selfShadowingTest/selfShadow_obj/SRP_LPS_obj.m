function sat = SRP_LPS_obj(sat, sun_b, d)
% ----------------------------------------------------------------------
%   calculate Solar Radiation Pressure (SRP) forces and torques
%                    with Lambertian and Perfect Specular model
%    20180816  y.yoshimura
%    Inputs: sat,
%             sun, sun vector from sat to sun, 3x1 unit vector
%               d, distance between sat and sun, m
%   Outputs: faceの中にforceとtorqueを加える
%   related function files:
%   note:
%   cf:
%   revisions;
%
%   (c) 2019 yasuhiro yoshimura
%----------------------------------------------------------------------
S0 = 1357; % Solar constant, W/m^2
c = 299792458; % light speed, m/s
AU = AU2km(1.0) * 10^3; % m
Bf = 2/3;
kappa = 0;

d_AU = d / AU; % AU
sun_b = sun_b(:); % column vector
sun_b = sun_b ./ norm(sun_b);

coeff = -S0 / c / d_AU^2;

sunlit_flag = (sat.normal * sun_b > 0); % 1: sunlit, 0: shade

tmp = abs(sat.normal*sun_b) .* (sat.Ca + sat.Cd) .* sun_b' ...
    + (sat.normal*sun_b) .* (Bf * sat.Cd + kappa * sat.Ca ...
    + 2.0 * sat.Cs * sat.normal * sun_b) .* sat.normal; % nx3 matrix
sat.force = sunlit_flag .* coeff .* sat.area .* tmp;
sat.torque = cross(sat.pos, sat.force);

end
