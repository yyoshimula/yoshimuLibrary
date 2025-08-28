%[text] # Calculate modified equinotical orbital element derivatives using Variational Equations
%[text] $\\frac{\\mathrm{d}p}{\\mathrm{d}t} = \\frac{2p}{w} \\sqrt{\\frac{p}{\\mu}} a\_T\\\\\n\\frac{\\mathrm{d}f}{\\mathrm{d}t} = \\sqrt{\\frac{p}{\\mu}}\\left( a\_R \\sin{L} + \\frac{((w+1)\\cos{L}+f)}{w}a\_T  - \\frac{g(h\\sin{L}-k\\cos{L})}{w}a\_N \\right)\\\\\n\\frac{\\mathrm{d}g}{\\mathrm{d}t} = \\sqrt{\\frac{p}{\\mu}}\\left( -a\_R \\cos{L} + \\frac{((w+1)\\sin{L}+g)}{w}a\_T  + \\frac{f(h\\sin{L}-k\\cos{L})}{w}a\_N\\right)  \\\\\n\\frac{\\mathrm{d h}}{\\mathrm{d}t} = \\sqrt{\\frac{p}{\\mu}}\\frac{s^2 a\_N}{2w}\\cos{L}\\\\\n\\frac{\\mathrm{d} k}{\\mathrm{d}t} = \\sqrt{\\frac{p}{\\mu}}\\frac{s^2 a\_N}{2w}\\sin{L}\\\\\n\\frac{\\mathrm{d}L}{\\mathrm{d}t} = \\sqrt{\\mu p }\\left(\\frac{w}{p}\\right)^2 + \\sqrt{\\frac{p}{\\mu}}\\frac{(h\\sin{L}-k\\cos{L})}{w}a\_N$
function dOEdt = mee(oe, aRTN, mu)

% Extract parameters for clarity
p = oe.p_;
f = oe.f_;
g = oe.g_;
h = oe.h_;
k = oe.k_;
L = oe.L_;

w = 1 + f * cos(L) + g * sin(L);
s2 = 1 + h^2 + k^2;

% Assemble Gauss Variational Equations
dOEdt = [ 2 * p / w * sqrt(p / mu) * aRTN(2)
    sqrt(p / mu) * (aRTN(1) * sin(L) + ((w + 1) * cos(L) + f) * aRTN(2) / w - g * (h * sin(L) - k * cos(L)) * aRTN(3) / w)
    sqrt(p / mu) * (-aRTN(1) * cos(L) + ((w + 1) * sin(L) + g)* aRTN(2) / w + f * (h * sin(L) - k * cos(L)) * aRTN(3) / w)
    sqrt(p / mu) * s2 * aRTN(3) * cos(L) / 2 / w
    sqrt(p / mu) * s2 * aRTN(3) * sin(L) / 2 / w
    sqrt(mu * p) * (w / p)^2 + sqrt(p / mu) * (h * sin(L) - k * cos(L)) * aRTN(3) / w
    ];
end

%[appendix]{"version":"1.0"}
%---
