%[text] # calculating light curves
%[text] light curves for the object defined with sat
%[text] `sat`: satellite configuration read with `readSC`
%[text] `scalar:` specify the definition of the quaternion 
%[text] `    scalar == 0`
%[text]         ${\\bf q} = \[q\_0,q\_1,q\_2,q\_3\]^T=\[\\cos(\\theta/2), {\\bf e}^T\\sin(\\theta/2)\]^T$
%[text] `    scalar == 4`
%[text]         ${\\bf q} = \[q\_1,q\_2,q\_3, q\_4\]^T=\[{\\bf e}^T\\sin(\\theta/2), \\cos(\\theta/2)\]^T$
%[text] `q`: quaternions, nx4 vector 
%[text] `satPos:`  satellite position vector@inertial frame, m, nx3 matrix
%[text] `obsPos: `observer position vector@inertial frame, m, nx3 matrix
%[text] `sunPos: sun` position vector@inertial frame, m, nx3 matrix
%[text] `nu:` flag if object is sunlit, nx1 vector
%[text]     `nu == 1` sunlit, `nu == 0` eclipse
%[text] BRDF, BRDF used in light curve clculation:
%[text]     default = Lambertian diffuse and mirror-like specular
%[text]        `AS` = Ashikhmin–Shirley model
%[text] `   CT` = Cook–Torrance model
%[text] `m`: relative magnitude of light curves w.r.t. Sun
%[text] $f\_{\\rm obs}\n$
%[text] ## note
%[text] 
%[text] ## references 
%[text] 
%[text] ## revisions
%[text] 20200430  y.yoshimura, y.yoshimula@gmail.com
%[text] See also readSC, lcSimple, lcAS, lcCT, orbitConst.
function [m, fObs] = lc(sat, scalar, q, satPos, obsPos, sunPos, nu, options)
arguments
    sat
    scalar (1,1) {mustBeMember(scalar, [0, 4])}
    q
    satPos
    obsPos
    sunPos
    nu
    options.BRDF (1,1) string {mustBeMember(options.BRDF, ["simple", "CT", "AS"])} = "simple"
    options.mex (1,1) string {mustBeMember(options.mex, ["on", "off"])} = "off"
end

%[text] ## Sun and observer relative position
sunRelDir = sunPos - satPos; % m, relative direction from sat to sun@ECI
obsRelDir = obsPos - satPos; % m, relative from sat to observer;

sunB = qRotation(scalar, normRow(sunRelDir), q); % sun direction (unit) vector@body-fixed frame, Nx3 mat
obsB = qRotation(scalar, normRow(obsRelDir), q); % observer direction (unit) vector@body-fixed frxame, Nx3 mat

%% relative magnitude of light curves
if strcmp(options.BRDF, 'simple')
    if strcmp(options.mex, 'on')
        sat = lcSimple_mex(sat, sunB, obsB);
    else
        sat = lcSimple(sat, sunB, obsB);
    end

elseif strcmp(options.BRDF, 'AS')
    if strcmp(options.mex, 'on')
        [sat, ~, ~] = lcAS_mex(sat, sunB, obsB);
    else
        [sat, ~, ~] = lcAS(sat, sunB, obsB);
    end
elseif strcmp(options.BRDF, 'CT') % Cook–Torrance model
    if strcmp(options.mex, 'on')
        sat = lcCT_mex(sat, sunB, obsB);
    else
        sat = lcCT(sat, sunB, obsB);
    end
else
    error('unavailable BRDF model')
end

fObs = sum(sat.fObs,1);
fObs = fObs(:); % column vector as time history
obs2Sun = normRow(sunPos - obsPos); % relative directional vector from obs to Sun
obs2Sat = normRow(satPos - obsPos); % relative directional vector from obs to sat

fObs = fObs .* nu; % consider umbra, penumbra
% fObs = fObs .* (dot(obs2Sun, normRow(obsPos),2) <= 0.0) % when observer cannot see Sun
% fObs = fObs .* (dot(obs2Sat, normRow(obsPos),2) > 0.0); % when observer can see Satellite

m = -26.7 - 2.5 * log10(fObs ./ vecnorm(obsRelDir,2,2).^2);

end

%[appendix]{"version":"1.0"}
%---
