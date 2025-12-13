%[text] Special case for Cook-Torrance SRP
%[text] $\\bf{s}=\\bf{n}$ where $\\bf n\n$ is along z-axis
clc
clear
cls

syms phi_r theta_r real
syms F0 m nRef real


sunB = [0; 0; 1];
n = [0; 0; 1];

v = [cos(phi_r) * sin(theta_r)
    sin(phi_r) * sin(theta_r)
    cos(theta_r)];

thetaH = theta_r / 2;
h = [cos(phi_r) * sin(thetaH)
    sin(phi_r) * sin(thetaH)
    cos(thetaH)];

% VH = v' * h;
VH = cos(thetaH);
NH = n' * h %[output:31a18e1a]
NV = n' * v %[output:5bf5a8ba]

g = sqrt(nRef^2 + VH^2 - 1);

D = exp(-tan(thetaH) / m)^2 / (pi * m^2 * cos(thetaH)^4);
% G = min(1, 2*NH/VH);
G = 2 * NH / VH;

tmp1 = (g - VH).^2 / 2 ./ (g + VH).^2;
tmp2 = (1 + (VH .* (g + VH) - 1).^2 ./ (VH .* (g - VH) + 1).^2);
F = tmp1 .* tmp2; % nFacet x N
F = F0;

cs = D * G * F / 4 / 1 / NV;

fx = cs * NV * sin(theta_r) .* v %[output:45987ffc]

simplify(int(int(fx(1), phi_r, 0, 2*pi), theta_r, 0, pi/2)) %[output:4614e60b]
simplify(int(int(fx(2), phi_r, 0, 2*pi), theta_r, 0, pi/2)) %[output:15edee67]
% simplify(int(int(fx(3), phi_r, 0, 2*pi), theta_r, 0, pi/2))
int(fx(3), theta_r, 0, pi/2) %[output:638aeb21]


%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
%[output:31a18e1a]
%   data: {"dataType":"symbolic","outputData":{"name":"NH","value":"\\cos \\left(\\frac{\\theta_r }{2}\\right)"}}
%---
%[output:5bf5a8ba]
%   data: {"dataType":"symbolic","outputData":{"name":"NV","value":"\\cos \\left(\\theta_r \\right)"}}
%---
%[output:45987ffc]
%   data: {"dataType":"symbolic","outputData":{"name":"fx","value":"\\begin{array}{l}\n\\left(\\begin{array}{c}\n\\frac{F_0 \\,{\\mathrm{e}}^{-\\frac{2\\,\\tan \\left(\\frac{\\theta_r }{2}\\right)}{m}} \\,\\cos \\left(\\phi_r \\right)\\,{\\sin \\left(\\theta_r \\right)}^2 }{\\sigma_1 }\\\\\n\\frac{F_0 \\,{\\mathrm{e}}^{-\\frac{2\\,\\tan \\left(\\frac{\\theta_r }{2}\\right)}{m}} \\,\\sin \\left(\\phi_r \\right)\\,{\\sin \\left(\\theta_r \\right)}^2 }{\\sigma_1 }\\\\\n\\frac{F_0 \\,{\\mathrm{e}}^{-\\frac{2\\,\\tan \\left(\\frac{\\theta_r }{2}\\right)}{m}} \\,\\cos \\left(\\theta_r \\right)\\,\\sin \\left(\\theta_r \\right)}{\\sigma_1 }\n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =2\\,m^2 \\,\\pi \\,{\\cos \\left(\\frac{\\theta_r }{2}\\right)}^4 \n\\end{array}"}}
%---
%[output:4614e60b]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"0"}}
%---
%[output:15edee67]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"0"}}
%---
%[output:638aeb21]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\int_0^{\\frac{\\pi }{2}} \\frac{F_0 \\,{\\mathrm{e}}^{-\\frac{2\\,\\tan \\left(\\frac{\\theta_r }{2}\\right)}{m}} \\,\\cos \\left(\\theta_r \\right)\\,\\sin \\left(\\theta_r \\right)}{2\\,m^2 \\,\\pi \\,{\\cos \\left(\\frac{\\theta_r }{2}\\right)}^4 }\\;\\textrm{d}\\theta_r"}}
%---
