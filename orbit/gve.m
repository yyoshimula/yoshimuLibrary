function dOEdt = gve(oe, aRTN, anomalyFlag)
    % Calculate orbital element derivatives using Gauss Variational Equations
    
    % Extract parameters for clarity
    r = oe.r;
    u = oe.u;
    n = oe.n;
    a = oe.a;
    e = oe.e;
    inc = oe.inc;
    f = oe.f;
    p = oe.p;
    h = oe.h;
    
    % Calculate RAAN rate
    dOmedt = r * sin(u) / n / a^2 / sqrt(1 - e^2) / sin(inc) * aRTN(3);
    
    % Calculate argument of perigee rate
    domedt = -sqrt(1 - e^2) / n / a / e * (cos(f) * aRTN(1) - ...
        (sin(f) + sin(f) / (1 + e * cos(f))) * aRTN(2)) - dOmedt * cos(inc);
    
    % Calculate anomaly rate
    if anomalyFlag == 1    % True anomaly
        tmp = h / r^2 - domedt - dOmedt * cos(inc);
    else                   % Mean anomaly
        tmp = n + 1 / n / a^2 / e * ((p * cos(f) - 2 * e * r) * aRTN(1) - ...
            (p + r) * sin(f) * aRTN(2));
    end
    
    % Assemble Gauss Variational Equations
    dOEdt = [
        2 / n / sqrt(1 - e^2) * (e * sin(f) * aRTN(1) + (1 + e * cos(f)) * aRTN(2))
        sqrt(1 - e^2) / n / a * (sin(f) * aRTN(1) + (cos(f) + (e + cos(f)) / (1 + e * cos(f))) * aRTN(2))
        r * cos(u) / (n * a^2 * sqrt(1 - e^2)) * aRTN(3)
        dOmedt
        domedt
        tmp
        ];
    end