%[text] # Calculate classical orbital element derivatives using Gauss Variational Equations
%[text] $\\oe = \[a, e, i, \\Omega, w, \\nu\]^T$ 
%[text] $\\frac{\\mathrm{d}a}{\\mathrm{d}t} = \\frac{2}{n\\sqrt{1 - e^{2}}}\\left\[e\\sin{\\nu}a\_R + (1+e\\cos{\\nu})a\_{T}\\right\] \\\\\n\\frac{\\mathrm{d}e}{\\mathrm{d}t} = \\frac{\\sqrt{1-e^{2}}}{na}\\left\[a\_R\\sin{\\nu} + \\left(\\cos{\\nu}+\\frac{e+\\cos{\\nu}}{1+e\\cos{\\nu}}\\right)a\_{T}\\right\] \\\\\n\\frac{\\mathrm{d}i}{\\mathrm{d}t} = \\frac{r \\cos{u}}{na^{2}\\sqrt{1-e^{2}}}a\_{N} \\\\\n\\frac{\\mathrm{d}\\Omega}{\\mathrm{d}t} = \\frac{r \\sin{u}}{na^{2}\\sqrt{1-e^{2}}\\sin{i}} a\_{N}\\\\\n\\frac{\\mathrm{d}w}{\\mathrm{d}t} = -\\frac{\\sqrt{1-e^{2}}}{nae}\\left\[a\_R\\cos{\\nu}-\\left(\\sin{\\nu}+\\frac{\\sin{\\nu}}{1+e\\cos{\\nu}}\\right)a\_{T}\\right\] - \\frac{\\mathrm{d}\\Omega}{\\mathrm{d}t}\\cos{i} \\\\\n\\frac{\\mathrm{d}\\nu}{\\mathrm{d}t} = \\frac{h}{r^{2}} -\\frac{\\mathrm{d}w}{\\mathrm{d}t}-\\frac{\\mathrm{d}\\Omega}{\\mathrm{d}t}\\cos{i}$
function dOEdt = gve(oe, aRTN, anomalyFlag, mu)
    oe = calcOrbitalState(oe, mu);
    r = oe.r;
    u = oe.u;
    n = oe.n;
    p = oe.p;
    h = oe.h;

    a = oe.a;
    e = oe.e;
    inc = oe.inc;
    nu = oe.nu;
    
    % Calculate RAAN rate
    dOmedt = r * sin(u) / n / a^2 / sqrt(1 - e^2) / sin(inc) * aRTN(3);
    
    % Calculate argument of perigee rate
    domedt = -sqrt(1 - e^2) / n / a / e * (cos(nu) * aRTN(1) - ...
        (sin(nu) + sin(nu) / (1 + e * cos(nu))) * aRTN(2)) - dOmedt * cos(inc);
    
    % Calculate anomaly rate
    if anomalyFlag == 1    % True anomaly
        tmp = h / r^2 - domedt - dOmedt * cos(inc);
    else                   % Mean anomaly assuming 
        tmp = n + 1 / n / a^2 / e * ((p * cos(nu) - 2 * e * r) * aRTN(1) - ...
            (p + r) * sin(nu) * aRTN(2));
    end
    
    % Assemble Gauss Variational Equations
    dOEdt = [
        2 / n / sqrt(1 - e^2) * (e * sin(nu) * aRTN(1) + (1 + e * cos(nu)) * aRTN(2))
        sqrt(1 - e^2) / n / a * (sin(nu) * aRTN(1) + (cos(nu) + (e + cos(nu)) / (1 + e * cos(nu))) * aRTN(2))
        r * cos(u) / (n * a^2 * sqrt(1 - e^2)) * aRTN(3)
        dOmedt
        domedt
        tmp
        ];
    end

%[appendix]{"version":"1.0"}
%---
