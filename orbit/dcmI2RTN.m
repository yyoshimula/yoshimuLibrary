%[text] # DCM from inertial frame to RTN frame
% Input: raan, inc, w, nu
% Output: R

function R = dcmI2RTN(raan, inc, w, nu)

R = zxz2dcm(raan, inc, w + nu);

end

%[appendix]{"version":"1.0"}
%---
