%[text] # Euler angle kinematics generator
%[text] ## parameters to be set
%[text] first: 1st rotation axis, must be 1, 2, or 3
%[text] second: 2nd rotation axis, must be 1, 2, or 3
%[text] third: 3rd rotation axis, must be 1, 2, or 3
%[text] ## output
%[text] dxdt: euler angle kinematics
%[text] ## note
%[text] NA
%[text] ## references 
%[text] Markley, F. L., and Crassidis, J. L., Fundamentals of Spacecraft Attitude Determination and Control, New York, NY: Springer, 2014. 
%[text] ## revisions
%[text] 20231212  y.yoshimura, y.yoshimula@gmail.com
%[text] See also .
clc
clear

%[text] ## to be set
first = 1;
second = 2;
third = 3;

%[text] ## derivation of matrix $B \\in \\mathbb{R}^{3\\times 3}$
%[text] $\\phi \\rightarrow \\theta \\rightarrow \\psi$の回転
%[text] $\\left\[\n\\begin{array}{c}\n    \\dot{\\phi}  \\\\\n     \\dot{\\theta}  \\\\\n     \\dot{\\psi}\n\\end{array}\\right\] = B(\\theta,\\psi)\\left\[\n\\begin{array}{c}\n    \\omega\_x  \\\\\n     \\omega\_y  \\\\\n     \\omega\_z\n\\end{array}\\right\]$
ePhi = zeros(3,1);
ePhi(first,1) = 1;

eTheta = zeros(3,1);
eTheta(second,1) = 1;

ePsi = zeros(3,1);
ePsi(third,1) = 1;

syms theta0
tmp = solve(ePsi == dcm1axis(second, theta0)*ePhi, theta0) %[output:1396486f]
theta0 = tmp;

syms theta psi real
B = [cross(ePsi, eTheta)'
    sin(theta - theta0) * eTheta'
    sin(theta - theta0) * ePsi' - cos(theta - theta0) * cross(ePsi, eTheta)'] * dcm1axis(third, psi)';
B = B ./ sin(theta - theta0);

B = simplify(B) %[output:88e8240a]

latex(B) %[output:26cde51b]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
%[output:1396486f]
%   data: {"dataType":"symbolic","outputData":{"name":"tmp","value":"\\frac{\\pi }{2}"}}
%---
%[output:88e8240a]
%   data: {"dataType":"symbolic","outputData":{"name":"B","value":"\\left(\\begin{array}{ccc}\n\\frac{\\cos \\left(\\psi \\right)}{\\cos \\left(\\theta \\right)} & -\\frac{\\sin \\left(\\psi \\right)}{\\cos \\left(\\theta \\right)} & 0\\\\\n\\sin \\left(\\psi \\right) & \\cos \\left(\\psi \\right) & 0\\\\\n-\\frac{\\cos \\left(\\psi \\right)\\,\\sin \\left(\\theta \\right)}{\\cos \\left(\\theta \\right)} & \\frac{\\sin \\left(\\psi \\right)\\,\\sin \\left(\\theta \\right)}{\\cos \\left(\\theta \\right)} & 1\n\\end{array}\\right)"}}
%---
%[output:26cde51b]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"'\\left(\\begin{array}{ccc} \\frac{\\cos\\left(\\psi \\right)}{\\cos\\left(\\theta \\right)} & -\\frac{\\sin\\left(\\psi \\right)}{\\cos\\left(\\theta \\right)} & 0\\\\ \\sin\\left(\\psi \\right) & \\cos\\left(\\psi \\right) & 0\\\\ -\\frac{\\cos\\left(\\psi \\right)\\,\\sin\\left(\\theta \\right)}{\\cos\\left(\\theta \\right)} & \\frac{\\sin\\left(\\psi \\right)\\,\\sin\\left(\\theta \\right)}{\\cos\\left(\\theta \\right)} & 1 \\end{array}\\right)'"}}
%---
