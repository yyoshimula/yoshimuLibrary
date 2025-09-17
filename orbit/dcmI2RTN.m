%[text] # DCM from inertial frame to RTN frame
% Input: raan, inc, w, f
% Output: R

function R = dcmI2RTN(raan, inc, w, f)

R = zxz2dcm(raan, inc, w + f);

end

%[appendix]{"version":"1.0"}
%---
