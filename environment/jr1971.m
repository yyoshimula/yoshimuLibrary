function result = jr1971(jd, phi_gd, lambda, h, F10, F10a, Kp)
% JR1971 Jacchia-Roberts 1971 Atmospheric Model
%
% Compute the atmospheric density using the Jacchia-Roberts 1971 model.
%
% Inputs:
%   jd     - Julian day
%   phi_gd - Geodetic latitude [rad]
%   lambda - Longitude [rad]
%   h      - Altitude [m]
%   F10    - 10.7-cm solar flux [sfu]
%   F10a   - 10.7-cm averaged solar flux, 81-day centered [sfu]
%   Kp     - Kp geomagnetic index (3-hour delayed)
%
%   F10, F10a, Kp は CelesTrak の SW-All.csv から lookupSolarGeoIndex() を
%   使って取得できる:
%       [F10, F10a, Kp] = lookupSolarGeoIndex(jd, 'SW-All.csv');
%   SW-All.csv は https://celestrak.org/SpaceData/SW-All.csv から入手可能。
%
% Output:
%   result - struct with fields:
%     total_density          [kg/m^3]
%     temperature            [K]
%     exospheric_temperature [K]
%     N2_number_density      [1/m^3]
%     O2_number_density      [1/m^3]
%     O_number_density       [1/m^3]
%     Ar_number_density      [1/m^3]
%     He_number_density      [1/m^3]
%     H_number_density       [1/m^3]
%
% Example:
%   % ISS 高度 (~400 km) における大気密度を計算する例
%   jd     = gc2jd(2024, 1, 1, 12, 0, 0); % 2024-01-01 12:00:00 UTC
%   phi_gd = deg2rad(35.0);   % 測地緯度 35 deg N
%   lambda = deg2rad(139.0);  % 経度 139 deg E
%   h      = 400e3;           % 高度 400 km → [m] で渡す
%   % F10/F10a/Kp は CelesTrak SW-All.csv から取得 (lookupSolarGeoIndex)
%   [F10, F10a, Kp] = lookupSolarGeoIndex(jd, 'SW-All.csv');
%   result = jr1971(jd, phi_gd, lambda, h, F10, F10a, Kp);
%   fprintf('rho = %.3e kg/m^3\n', result.total_density);
%
% References:
%   [1] Roberts, C. R (1971). An analytic model for upper atmosphere
%       densities based upon Jacchia's 1970 models.
%   [2] Jacchia, L. G (1970). New static models of the thermosphere and
%       exosphere with empirical temperature profiles. SAO Special Report #313.
%   [3] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
%       4th ed. Microcosm Press.
%   [4] Long, A. C. et al. (1989). GTDS Mathematical Theory (Revision 1).

    % == Constants =============================================================
    persistent C_cached
    if isempty(C_cached), C_cached = jr1971_constants(); end
    C = C_cached;

    Rstar = C.Rstar;
    Av    = C.Av;
    Ra    = C.Ra;
    g0    = C.g0;
    M0    = C.M0;
    z1    = C.z1;
    z2    = C.z2;
    T1    = C.T1;
    M1    = C.M1;
    rho1  = C.rho1;
    zx    = C.zx;
    Mi    = C.Mi;
    alphai = C.alphai;
    mui   = C.mui;
    Aa    = C.Aa;
    Ca    = C.Ca;
    la    = C.la;
    alpha = C.alpha;
    beta  = C.beta;
    zeta  = C.zeta;
    deltaij = C.deltaij;

    % == Auxiliary variables ===================================================
    Ra2 = Ra * Ra;

    % == Preliminaries =========================================================

    % Convert altitude from [m] to [km].
    h = h / 1000;

    % Compute the Sun position in MOD frame.
    s_i = sun_position_mod(jd);

    % Sun declination [rad].
    delta_s = atan2(s_i(3), sqrt(s_i(1)^2 + s_i(2)^2));

    % Sun right ascension [rad].
    Omega_s = atan2(s_i(2), s_i(1));

    % Right ascension of the selected location w.r.t. inertial frame.
    Omega_p = lambda + jd_to_gmst(jd);

    % Hour angle.
    H = Omega_p - Omega_s;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                               Algorithm                                  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % == Exospheric Temperature ================================================

    % -- Diurnal Variation -----------------------------------------------------

    dF10 = F10 - F10a;
    Tc   = 379 + 3.24 * F10 + 1.3 * dF10;

    eta   = abs(phi_gd - delta_s) / 2;
    theta = abs(phi_gd + delta_s) / 2;

    tau = wrapToPi(H + deg2rad(-37 + 6 * sin(H + deg2rad(43))));

    Cv = cos(eta)^2.2;
    S  = sin(theta)^2.2;
    Tl = Tc * (1 + 0.3 * (S + (Cv - S) * cos(tau / 2)^3));

    % == Geomagnetic Activity ==================================================

    if h < 200
        dTinf = 14 * Kp + 0.02 * exp(Kp);
    else
        dTinf = 28 * Kp + 0.03 * exp(Kp);
    end

    Tinf = Tl + dTinf;

    % == Temperature at the Desired Altitude ===================================

    a_coeff =  371.6678;
    b_coeff =    0.0518806;
    c_coeff = -294.3505;
    d_coeff =   -0.00216222;
    Tx = a_coeff + b_coeff * Tinf + c_coeff * exp(d_coeff * Tinf);

    Tz = jr1971_temperature(h, Tx, Tinf, C);

    % == Corrections to the Density ============================================

    % -- Geomagnetic Effect ----------------------------------------------------
    if h < 200
        dlog10rho_g = 0.012 * Kp + 1.2e-5 * exp(Kp);
    else
        dlog10rho_g = 0.0;
    end

    % -- Semi-annual Variation -------------------------------------------------
    Phi = (jd - 2436204.5) / 365.2422;

    tau_sa = Phi + 0.09544 * ((0.5 * (1 + sin(2*pi*Phi + 6.035)))^1.65 - 0.5);
    f_z    = (5.876e-7 * h^2.331 + 0.06328) * exp(-0.002868 * h);
    g_t    = 0.02835 + (0.3817 + 0.17829 * sin(2*pi*tau_sa + 4.137)) ...
             * sin(4*pi*tau_sa + 4.259);

    dlog10rho_sa = f_z * g_t;

    % -- Seasonal Latitudinal Variation ----------------------------------------
    sin_phi_gd     = sin(phi_gd);
    abs_sin_phi_gd = abs(sin_phi_gd);

    dlog10rho_lt = 0.014 * (h - 90) * exp(-0.0013 * (h - 90)^2) ...
                   * sin(2*pi*Phi + 1.72) * sin_phi_gd * abs_sin_phi_gd;

    % -- Total Correction ------------------------------------------------------
    dlog10rho_c = dlog10rho_g + dlog10rho_lt + dlog10rho_sa;
    drho_c      = 10^dlog10rho_c;

    % == Density ===============================================================

    if h == z1
        rho = rho1 * drho_c;

        result = make_output( ...
            1000*rho, Tz, Tinf, ...
            (rho * mui(1)) * Av / Mi(1) * 1e6, ...
            (rho * mui(2)) * Av / Mi(2) * 1e6, ...
            (rho * mui(3)) * Av / Mi(3) * 1e6, ...
            (rho * mui(4)) * Av / Mi(4) * 1e6, ...
            (rho * mui(5)) * Av / Mi(5) * 1e6, ...
            (rho * mui(6)) * Av / Mi(6) * 1e6);

    elseif z1 < h && h <= zx

        % Roots of the polynomial P(Z).
        c0 = (35^4 * Tx / (Tx - T1) + Ca(1)) / Ca(5);
        c1 = Ca(2) / Ca(5);
        c2 = Ca(3) / Ca(5);
        c3 = Ca(4) / Ca(5);
        c4 = 1.0;  % Ca(5)/Ca(5)

        [r1, r2, xr, yr] = jr1971_roots([c0, c1, c2, c3, c4]);

        % -- f and k -------------------------------------------------------
        f = 35^4 * Ra2 / Ca(5);
        k = -g0 / (Rstar * (Tx - T1));

        % -- X -------------------------------------------------------------
        X = -2 * r1 * r2 * Ra * (Ra2 + 2*xr*Ra + xr^2 + yr^2);

        if h <= z2

            % == Altitudes Between 90 km and 100 km ========================

            B0 = alpha(1) + beta(1) * Tx / (Tx - T1);
            B1 = alpha(2) + beta(2) * Tx / (Tx - T1);
            B2 = alpha(3) + beta(3) * Tx / (Tx - T1);
            B3 = alpha(4) + beta(4) * Tx / (Tx - T1);
            B4 = alpha(5) + beta(5) * Tx / (Tx - T1);
            B5 = alpha(6) + beta(6) * Tx / (Tx - T1);

            p2 =  evalpoly_asc(r1, [B0,B1,B2,B3,B4,B5]) / U_func(r1, Ra, xr, yr, r1, r2);
            p3 = -evalpoly_asc(r2, [B0,B1,B2,B3,B4,B5]) / U_func(r2, Ra, xr, yr, r1, r2);
            p5 =  evalpoly_asc(-Ra, [B0,B1,B2,B3,B4,B5]) / V_func(-Ra, xr, yr, r1, r2);

            p4 = (B0 - r1*r2*Ra2*(B4 + (2*xr + r1 + r2 - Ra)*B5) ...
                  - r1*r2*Ra*(xr^2 + yr^2)*B5 ...
                  + r1*r2*(Ra2 - (xr^2 + yr^2))*p5 ...
                  + W_func(r1, Ra, xr, yr, r1, r2)*p2 ...
                  + W_func(r2, Ra, xr, yr, r1, r2)*p3) / X;

            p6 = B4 + (2*xr + r1 + r2 - Ra)*B5 - p5 ...
                 - 2*(xr + Ra)*p4 - (r2 + Ra)*p3 - (r1 + Ra)*p2;
            p1 = B5 - 2*p4 - p3 - p2;

            % -- F1 and F2 -------------------------------------------------
            log_F1 = p1 * log((h + Ra) / (z1 + Ra)) ...
                   + p2 * log((h - r1) / (z1 - r1)) ...
                   + p3 * log((h - r2) / (z1 - r2)) ...
                   + p4 * log((h^2 - 2*xr*h + xr^2 + yr^2) / ...
                              (z1^2 - 2*xr*z1 + xr^2 + yr^2));

            F2 = (h - z1) * (Aa(7) + p5 / ((h + Ra) * (z1 + Ra))) ...
               + p6 / yr * atan(yr * (h - z1) / ...
                                (yr^2 + (h - xr) * (z1 - xr)));

            % -- Density ---------------------------------------------------
            Mz  = jr1971_mean_molecular_mass(h, Aa);
            rho = rho1 * drho_c * Mz * T1 / (M1 * Tz) * exp(k * (log_F1 + F2));

            result = make_output( ...
                1000*rho, Tz, Tinf, ...
                (rho * mui(1)) * Av / Mi(1) * 1e6, ...
                (rho * mui(2)) * Av / Mi(2) * 1e6, ...
                (rho * mui(3)) * Av / Mi(3) * 1e6, ...
                (rho * mui(4)) * Av / Mi(4) * 1e6, ...
                (rho * mui(5)) * Av / Mi(5) * 1e6, ...
                (rho * mui(6)) * Av / Mi(6) * 1e6);

        else

            % == Altitudes Between 100 km and 125 km =======================

            T100 = jr1971_temperature(z2, Tx, Tinf, C);

            rho100 = evalpoly_asc(Tinf, zeta) * M0;
            rho100 = rho100 * drho_c;

            q2 =  1 / U_func(r1, Ra, xr, yr, r1, r2);
            q3 = -1 / U_func(r2, Ra, xr, yr, r1, r2);
            q5 =  1 / V_func(-Ra, xr, yr, r1, r2);
            q4 = (1 + r1*r2*(Ra2 - (xr^2 + yr^2))*q5 ...
                  + W_func(r1, Ra, xr, yr, r1, r2)*q2 ...
                  + W_func(r2, Ra, xr, yr, r1, r2)*q3) / X;
            q6 = -q5 - 2*(xr + Ra)*q4 - (r2 + Ra)*q3 - (r1 + Ra)*q2;
            q1 = -2*q4 - q3 - q2;

            log_F3 = q1 * log((h + Ra) / (z2 + Ra)) ...
                   + q2 * log((h - r1) / (z2 - r1)) ...
                   + q3 * log((h - r2) / (z2 - r2)) ...
                   + q4 * log((h^2 - 2*xr*h + xr^2 + yr^2) / ...
                              (z2^2 - 2*xr*z2 + xr^2 + yr^2));

            F4 = q5 * (h - z2) / ((h + Ra) * (Ra + z2)) ...
               + q6 / yr * atan(yr * (h - z2) / ...
                                (yr^2 + (h - xr) * (z2 - xr)));

            expk = k * f * (log_F3 + F4);
            rhoN2 = rho100 * Mi(1)/M0 * mui(1) * (T100/Tz)^(1 + alphai(1)) * exp(Mi(1)*expk);
            rhoO2 = rho100 * Mi(2)/M0 * mui(2) * (T100/Tz)^(1 + alphai(2)) * exp(Mi(2)*expk);
            rhoO  = rho100 * Mi(3)/M0 * mui(3) * (T100/Tz)^(1 + alphai(3)) * exp(Mi(3)*expk);
            rhoAr = rho100 * Mi(4)/M0 * mui(4) * (T100/Tz)^(1 + alphai(4)) * exp(Mi(4)*expk);
            rhoHe = rho100 * Mi(5)/M0 * mui(5) * (T100/Tz)^(1 + alphai(5)) * exp(Mi(5)*expk);

            rho = rhoN2 + rhoO2 + rhoO + rhoAr + rhoHe;

            result = make_output( ...
                1000*rho, Tz, Tinf, ...
                rhoN2 * Av / Mi(1) * 1e6, ...
                rhoO2 * Av / Mi(2) * 1e6, ...
                rhoO  * Av / Mi(3) * 1e6, ...
                rhoAr * Av / Mi(4) * 1e6, ...
                rhoHe * Av / Mi(5) * 1e6, ...
                0.0);
        end

    else

        % == Altitudes Higher than 125 km ======================================

        rho125_N2 = drho_c * Mi(1) * 10^evalpoly_asc(Tinf, deltaij.N2) / Av;
        rho125_O2 = drho_c * Mi(2) * 10^evalpoly_asc(Tinf, deltaij.O2) / Av;
        rho125_O  = drho_c * Mi(3) * 10^evalpoly_asc(Tinf, deltaij.O)  / Av;
        rho125_Ar = drho_c * Mi(4) * 10^evalpoly_asc(Tinf, deltaij.Ar) / Av;
        rho125_He = drho_c * Mi(5) * 10^evalpoly_asc(Tinf, deltaij.He) / Av;

        % -- Compute l ---------------------------------------------------------
        l = evalpoly_asc(Tinf, la);

        % -- Eq. 25' -----------------------------------------------------------
        gamma   = (g0 * Ra2 / (Rstar * l * Tinf) * (Tinf - Tx) / (Tx - T1) ...
                   * (zx - z1) / (Ra + zx));
        gammaN2 = gamma * Mi(1);
        gammaO2 = gamma * Mi(2);
        gammaO  = gamma * Mi(3);
        gammaAr = gamma * Mi(4);
        gammaHe = gamma * Mi(5);

        % -- Eq. 25 ------------------------------------------------------------
        rhoN2 = rho125_N2 * (Tx/Tz)^(1 + alphai(1) + gammaN2) ...
                * ((Tinf - Tz)/(Tinf - Tx))^gammaN2;
        rhoO2 = rho125_O2 * (Tx/Tz)^(1 + alphai(2) + gammaO2) ...
                * ((Tinf - Tz)/(Tinf - Tx))^gammaO2;
        rhoO  = rho125_O  * (Tx/Tz)^(1 + alphai(3) + gammaO) ...
                * ((Tinf - Tz)/(Tinf - Tx))^gammaO;
        rhoAr = rho125_Ar * (Tx/Tz)^(1 + alphai(4) + gammaAr) ...
                * ((Tinf - Tz)/(Tinf - Tx))^gammaAr;
        rhoHe = rho125_He * (Tx/Tz)^(1 + alphai(5) + gammaHe) ...
                * ((Tinf - Tz)/(Tinf - Tx))^gammaHe;

        % -- Helium Seasonal/Latitudinal Correction ----------------------------
        dlog10rho_He = 0.65 / deg2rad(23.439291) * abs(delta_s) ...
                       * (sin(pi/4 - phi_gd*delta_s/(2*abs(delta_s)))^3 - 0.35355);
        rhoHe = rhoHe * 10^dlog10rho_He;

        % -- H for Altitude > 500 km -------------------------------------------
        rhoH = 0.0;
        if h > 500
            T500       = jr1971_temperature(500.0, Tx, Tinf, C);
            log10_T500 = log10(T500);
            rho500_H   = Mi(6) / Av * 10^(73.13 - (39.4 - 5.5*log10_T500)*log10_T500);

            gammaH = Mi(6) * gamma;
            rhoH   = drho_c * rho500_H * (T500/Tz)^(1 + gammaH) ...
                     * ((Tinf - Tz)/(Tinf - T500))^gammaH;
        end

        rho = rhoN2 + rhoO2 + rhoO + rhoAr + rhoHe + rhoH;

        result = make_output( ...
            1000*rho, Tz, Tinf, ...
            rhoN2 * Av / Mi(1) * 1e6, ...
            rhoO2 * Av / Mi(2) * 1e6, ...
            rhoO  * Av / Mi(3) * 1e6, ...
            rhoAr * Av / Mi(4) * 1e6, ...
            rhoHe * Av / Mi(5) * 1e6, ...
            rhoH  * Av / Mi(6) * 1e6);
    end
end

%% =========================================================================
%  Local Functions
%  =========================================================================

function result = make_output(total_density, temperature, exospheric_temperature, ...
                              N2, O2, O, Ar, He, H)
    result.total_density          = total_density;
    result.temperature            = temperature;
    result.exospheric_temperature = exospheric_temperature;
    result.N2_number_density      = N2;
    result.O2_number_density      = O2;
    result.O_number_density       = O;
    result.Ar_number_density      = Ar;
    result.He_number_density      = He;
    result.H_number_density       = H;
end

function C = jr1971_constants()
    C.Rstar = 8.31432;       % Universal gas constant [J/(K*mol)]
    C.Av    = 6.022045e23;   % Avogadro's constant [molecules/mol]
    C.Ra    = 6356.766;      % Mean Earth radius [km]
    C.g0    = 9.80665;       % Gravity at sea level [m/s^2]

    C.M0 = 28.960;          % Mean molecular mass at sea level [g/mol]

    C.T1   = 183.0;         % Temperature at lower bound [K]
    C.z1   = 90.0;          % Altitude of lower bound [km]
    C.M1   = 28.82678;      % Mean molecular mass at z1 [g/mol]
    C.rho1 = 3.46e-9;       % Density at z1 [g/cm^3]

    C.z2 = 100.0;           % Second division altitude [km]
    C.zx = 125.0;           % Inflection point altitude [km]

    % Molecular mass [g/mol]: N2, O2, O, Ar, He, H
    C.Mi = [28.0134, 31.9988, 15.9994, 39.9480, 4.0026, 1.00797];

    % Thermal diffusion coefficient: N2, O2, O, Ar, He, H
    C.alphai = [0, 0, 0, 0, -0.38, 0];

    % Fractional volume composition: N2, O2, O, Ar, He, H
    C.mui = [0.78110, 0.161778, 0.095544, 0.0093432, 0.61471e-5, 0];

    C.Aa = [-435093.363387, ...
              28275.5646391, ...
               -765.33466108, ...
                 11.043387545, ...
                 -0.08958790995, ...
                  0.00038737586, ...
                 -0.000000697444];

    C.Ca = [-89284375.0, ...
              3542400.0, ...
               -52687.5, ...
                  340.5, ...
                   -0.8];

    C.la = [0.1031445e+5, ...
            0.2341230e+1, ...
            0.1579202e-2, ...
           -0.1252487e-5, ...
            0.2462708e-9];

    C.alpha = [3144902516.672729, ...
              -123774885.4832917, ...
                 1816141.096520398, ...
                  -11403.31079489267, ...
                      24.36498612105595, ...
                       0.008957502869707995];

    C.beta = [-52864482.17910969, ...
                -16632.50847336828, ...
                    -1.308252378125, ...
                     0, 0, 0];

    C.zeta = [0.1985549e-10, ...
             -0.1833490e-14, ...
              0.1711735e-17, ...
             -0.1021474e-20, ...
              0.3727894e-24, ...
             -0.7734110e-28, ...
              0.7026942e-32];

    C.deltaij.N2 = [0.1093155e2, ...
                    0.1186783e-2, ...
                   -0.1677341e-5, ...
                    0.1420228e-8, ...
                   -0.7139785e-12, ...
                    0.1969715e-15, ...
                   -0.2296182e-19];

    C.deltaij.O2 = [0.9924237e1, ...
                    0.1600311e-2, ...
                   -0.2274761e-5, ...
                    0.1938454e-8, ...
                   -0.9782183e-12, ...
                    0.2698450e-15, ...
                   -0.3131808e-19];

    C.deltaij.O  = [0.1097083e2, ...
                    0.6118742e-4, ...
                   -0.1165003e-6, ...
                    0.9239354e-10, ...
                   -0.3490739e-13, ...
                    0.5116298e-17, ...
                    0];

    C.deltaij.Ar = [0.8049405e1, ...
                    0.2382822e-2, ...
                   -0.3391366e-5, ...
                    0.2909714e-8, ...
                   -0.1481702e-11, ...
                    0.4127600e-15, ...
                   -0.4837461e-19];

    C.deltaij.He = [0.7646886e1, ...
                   -0.4383486e-3, ...
                    0.4694319e-6, ...
                   -0.2894886e-9, ...
                    0.9451989e-13, ...
                   -0.1270838e-16, ...
                    0];
end

function T = jr1971_temperature(z, Tx, Tinf, C)
% Compute the temperature [K] at height z [km].
    T1 = C.T1;
    z1 = C.z1;
    zx = C.zx;

    if z < z1
        error('The altitude must not be lower than %.0f km.', z1);
    end
    if Tinf < 0
        error('The exospheric temperature must be positive.');
    end

    if z <= zx
        Ca  = C.Ca;
        aux = evalpoly_asc(z, Ca);
        T   = Tx + (Tx - T1) / 35^4 * aux;
    else
        Ra = C.Ra;
        la = C.la;
        l  = evalpoly_asc(Tinf, la);
        T  = Tinf - (Tinf - Tx) * exp( ...
             -l * ((Tx - T1)/(Tinf - Tx)) * ((z - zx)/(zx - z1)) / (Ra + z));
    end
end

function Mz = jr1971_mean_molecular_mass(z, Aa)
% Compute mean molecular mass at altitude z [km] (valid for 90 <= z <= 100 km).
    if z < 90 || z > 100
        warning('Empirical mean molecular mass model valid only for 90 <= z <= 100 km.');
    end
    Mz = evalpoly_asc(z, Aa);
end

function [r1, r2, xr, yr] = jr1971_roots(coeffs)
% Compute roots of the 4th-degree polynomial for density below 125 km.
% coeffs = [c0, c1, c2, c3, c4] (ascending order).
% MATLAB's roots() expects descending order.
    p_desc = fliplr(coeffs);
    r = roots(p_desc);

    % Separate real and complex roots.
    tol = 1e-10;
    real_mask = abs(imag(r)) < tol;
    real_roots = real(r(real_mask));
    complex_roots = r(~real_mask);

    r1 = max(real_roots);
    r2 = min(real_roots);

    % Find the complex root with positive imaginary part.
    idx = find(imag(complex_roots) > 0, 1);
    xr = real(complex_roots(idx));
    yr = abs(imag(complex_roots(idx)));
end

function val = U_func(nu, Ra, x, y, r1, r2)
    val = (nu + Ra)^2 * (nu^2 - 2*x*nu + x^2 + y^2) * (r1 - r2);
end

function val = V_func(nu, x, y, r1, r2)
    val = (nu^2 - 2*x*nu + x^2 + y^2) * (nu - r1) * (nu - r2);
end

function val = W_func(nu, Ra, x, y, r1, r2)
    val = r1 * r2 * Ra * (Ra + nu) * (Ra + (x^2 + y^2)/nu);
end

function val = evalpoly_asc(x, coeffs)
% Evaluate polynomial with coefficients in ascending order (c0, c1, c2, ...).
% p(x) = c0 + c1*x + c2*x^2 + ... (Horner's method)
    n = length(coeffs);
    val = coeffs(n);
    for i = n-1:-1:1
        val = val * x + coeffs(i);
    end
end

function s = sun_position_mod(jd)
% Compute the Sun position vector in the Mean-of-Date (MOD) frame [km].
% Low-precision solar coordinates (Vallado, Meeus).
    T_UT1 = (jd - 2451545.0) / 36525.0;

    % Mean longitude of the Sun [deg].
    lambda_M = mod(280.46 + 36000.771 * T_UT1, 360);

    % Mean anomaly of the Sun [deg].
    M_sun = mod(357.5291092 + 35999.0502909 * T_UT1, 360);
    M_rad = deg2rad(M_sun);

    % Ecliptic longitude [deg].
    lambda_ec = lambda_M + 1.914666471 * sin(M_rad) + 0.019994643 * sin(2*M_rad);
    lambda_ec_rad = deg2rad(lambda_ec);

    % Obliquity of the ecliptic [deg].
    epsilon = 23.439291 - 0.0130042 * T_UT1;
    epsilon_rad = deg2rad(epsilon);

    % Distance to the Sun [AU].
    r_sun = 1.000140612 - 0.016708617 * cos(M_rad) - 0.000139589 * cos(2*M_rad);

    % Convert to km (1 AU = 149597870.7 km).
    r_km = r_sun * 149597870.7;

    % Sun position in MOD frame [km].
    s = [r_km * cos(lambda_ec_rad); ...
         r_km * cos(epsilon_rad) * sin(lambda_ec_rad); ...
         r_km * sin(epsilon_rad) * sin(lambda_ec_rad)];
end

function gmst = jd_to_gmst(jd)
% Compute Greenwich Mean Sidereal Time [rad] from Julian Date.
    T_UT1 = (jd - 2451545.0) / 36525.0;

    % GMST in seconds.
    gmst_sec = 67310.54841 ...
             + (876600*3600 + 8640184.812866) * T_UT1 ...
             + 0.093104 * T_UT1^2 ...
             - 6.2e-6 * T_UT1^3;

    % Convert to radians (mod 2*pi).
    gmst = mod(gmst_sec * pi / 43200, 2*pi);
end
