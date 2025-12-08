clc
cls
clear

%[text] integration
syms theta_h phi_h theta_n lambda_c real

s = [0; 0; 1];

vh = [sin(2*theta_h) * cos(phi_h)
    sin(2*theta_h) * sin(phi_h)
    cos(2*theta_h)];
h = [sin(theta_h) * cos(phi_h) %[output:group:03c0cce3] %[output:2a925b98]
    sin(theta_h) * sin(phi_h) %[output:2a925b98]
    cos(theta_h)] %[output:group:03c0cce3] %[output:2a925b98]

n = [0
    sin(theta_n)
    cos(theta_n)];

sg1 = exp(lambda_c*(vh'*n - 1));
integrand = 2 * sg1 * sin(2 * theta_h) .* vh %[output:1b23c2dc]

A = int(int(integrand, phi_h), theta_h) %[output:6791cb2b]
simplify(A(2:3)) %[output:2b791992]

theta_h = 0.2;
phi_h = 0.1;

phi_h * sin(2*theta_h)^2 / 2 %[output:0017434b]

-1/4*phi_h * cos(4*theta_h) %[output:7b3f5d7d]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":55.8}
%---
%[output:2a925b98]
%   data: {"dataType":"symbolic","outputData":{"name":"h","value":"\\left(\\begin{array}{c}\n\\cos \\left(\\phi_h \\right)\\,\\sin \\left(\\theta_h \\right)\\\\\n\\sin \\left(\\phi_h \\right)\\,\\sin \\left(\\theta_h \\right)\\\\\n\\cos \\left(\\theta_h \\right)\n\\end{array}\\right)"}}
%---
%[output:1b23c2dc]
%   data: {"dataType":"symbolic","outputData":{"name":"integrand","value":"\\begin{array}{l}\n\\left(\\begin{array}{c}\n2\\,{\\sin \\left(2\\,\\theta_h \\right)}^2 \\,\\sigma_1 \\,\\cos \\left(\\phi_h \\right)\\\\\n2\\,{\\sin \\left(2\\,\\theta_h \\right)}^2 \\,\\sigma_1 \\,\\sin \\left(\\phi_h \\right)\\\\\n2\\,\\cos \\left(2\\,\\theta_h \\right)\\,\\sin \\left(2\\,\\theta_h \\right)\\,\\sigma_1 \n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 ={\\mathrm{e}}^{\\lambda_c \\,{\\left(\\cos \\left(2\\,\\theta_h \\right)\\,\\cos \\left(\\theta_n \\right)+\\sin \\left(2\\,\\theta_h \\right)\\,\\sin \\left(\\phi_h \\right)\\,\\sin \\left(\\theta_n \\right)-1\\right)}} \n\\end{array}"}}
%---
%[output:6791cb2b]
%   data: {"dataType":"symbolic","outputData":{"name":"A","value":"\\begin{array}{l}\n\\left(\\begin{array}{c}\n\\int \\frac{4\\,{\\mathrm{e}}^{-\\lambda_c \\,\\cos \\left(\\theta_n \\right)} \\,{\\mathrm{e}}^{-\\lambda_c } \\,{\\mathrm{e}}^{2\\,\\lambda_c \\,{\\cos \\left(\\theta_h \\right)}^2 \\,\\cos \\left(\\theta_n \\right)} \\,{\\mathrm{e}}^{2\\,\\lambda_c \\,\\cos \\left(\\theta_h \\right)\\,\\sin \\left(\\phi_h \\right)\\,\\sin \\left(\\theta_h \\right)\\,\\sin \\left(\\theta_n \\right)} \\,\\cos \\left(\\theta_h \\right)\\,\\sin \\left(\\theta_h \\right)}{\\lambda_c \\,\\sin \\left(\\theta_n \\right)}\\textrm{d}\\theta_h \\\\\n\\int \\int 2\\,{\\sin \\left(2\\,\\theta_h \\right)}^2 \\,\\sigma_1 \\,\\sin \\left(\\phi_h \\right)\\textrm{d}\\phi_h \\textrm{d}\\theta_h \\\\\n\\int \\int 2\\,\\cos \\left(2\\,\\theta_h \\right)\\,\\sin \\left(2\\,\\theta_h \\right)\\,\\sigma_1 \\textrm{d}\\phi_h \\textrm{d}\\theta_h \n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 ={\\mathrm{e}}^{\\lambda_c \\,{\\left(\\cos \\left(2\\,\\theta_h \\right)\\,\\cos \\left(\\theta_n \\right)+\\sin \\left(2\\,\\theta_h \\right)\\,\\sin \\left(\\phi_h \\right)\\,\\sin \\left(\\theta_n \\right)-1\\right)}} \n\\end{array}"}}
%---
%[output:2b791992]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\begin{array}{l}\n\\left(\\begin{array}{c}\n\\int \\int 2\\,{\\sin \\left(2\\,\\theta_h \\right)}^2 \\,\\sigma_1 \\,\\sin \\left(\\phi_h \\right)\\textrm{d}\\phi_h \\textrm{d}\\theta_h \\\\\n\\int \\int 2\\,\\cos \\left(2\\,\\theta_h \\right)\\,\\sin \\left(2\\,\\theta_h \\right)\\,\\sigma_1 \\textrm{d}\\phi_h \\textrm{d}\\theta_h \n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 ={\\mathrm{e}}^{\\lambda_c \\,{\\left(\\cos \\left(2\\,\\theta_h \\right)\\,\\cos \\left(\\theta_n \\right)+\\sin \\left(2\\,\\theta_h \\right)\\,\\sin \\left(\\phi_h \\right)\\,\\sin \\left(\\theta_n \\right)-1\\right)}} \n\\end{array}"}}
%---
%[output:0017434b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"0.0076"},"version":0}
%---
%[output:7b3f5d7d]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"-0.0174"},"version":0}
%---
