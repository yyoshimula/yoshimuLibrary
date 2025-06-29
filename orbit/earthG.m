%[text] # earth gravitational force w.r.t. ECEF frame
%[text] ## inputs
%[text] `jd`: Julian day, day, 
%[text] `rVec`: satellite position at inertial frame, km, nx3 vector
%[text] `const`: orbital constants
%[text] `EGM`: earth gravity constants
%[text] `options`
%[text]  `SPICE`: use SPICE, `on` or `off`
%[text] ## output
%[text] `aEarth`: Earth's gravitational force w.r.t. ECEF frame, km/s^2, nx3 vector
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20211027  y.yoshimura, y.yoshimula@gmail.com
%[text] See also orbitConst, vsopConst, precession.
function aEarth = earthG(jd, rVec, const, EGM, options)
arguments
    jd (:,1) {mustBeNumeric}
    rVec (:,3) {mustBeNumeric}
    const   
    EGM
    options.SPICE char {mustBeMember(options.SPICE,{'on', 'off'})} = 'off'
end

%[text] ### default: SPICE is not used
if strcmp(options.SPICE, 'off')
    dcmNutation = nutationDCM(jd,const);% nutation
    dcmPrecession = precessionDCM(const.J2000, jd, const); % precession
    GAST = gast(jd, const);
    i2PEF = dcm1axis(3, GAST) * dcmNutation * dcmPrecession;

    rPEF = i2PEF * rVec';
    
    tmp = egm2008(rPEF, EGM.GEODEG, EGM.Cnm, EGM.Snm, const); % at Cartesian coordinate @ECEF frame

%[text] ### when SPICE used
else
    % Julian day to ephemeris time (et)
    et = cspice_unitim(jd, 'JDTDT', 'ET');

    % state transition matrix from J2000 to ITRF93(ECEF)
    stm = cspice_sxform('J2000', 'ITRF93', et); % 6x6 matrix

    rITRF = stm(1:3,1:3) * rVec';

    tmp = egm2008(rITRF, EGM.GEODEG, EGM.Cnm, EGM.Snm, const); % at Cartesian coordinate

end

aEarth = tmp;

end

%[appendix]{"version":"1.0"}
%---
