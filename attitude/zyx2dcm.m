%[text] # ZYX Euler angle to directional cosine matrix (rotation matrix)
%[text] ## inputs
%[text] `phi`: rotation angle around Z-axis (1st rotation)
%[text] `theta`: rotation angle around Y-axis (2nd rotation)
%[text] `psi`: rotation angle around X-axis (3rd rotation)
%[text] ## output
%[text] R: directional cosine matrix (rotation matrix), 3x3
%[text] ## note
%[text] NA
%[text] ## references 
%[text] Markley, F. L., and Crassidis, J. L., Fundamentals of Spacecraft Attitude Determination and Control, New York, NY: Springer, 2014. 
%[text] ## revisions
%[text] 20150101  y.yoshimura, y.yoshimula@gmail.com
%[text] See also q2DCM.
function R = zyx2dcm(phi, theta, psi)
% arguments
%     phi (1,1)
%     theta (1,1)
%     psi (1,1) 
% end

R = [cos(theta)*cos(phi) cos(theta)*sin(phi) -sin(theta)
    sin(theta)*cos(phi)*sin(psi)-sin(phi)*cos(psi) sin(theta)*sin(phi)*sin(psi)+cos(phi)*cos(psi) cos(theta)*sin(psi)
    sin(theta)*cos(phi)*cos(psi)+sin(phi)*sin(psi) sin(theta)*sin(phi)*cos(psi)-cos(phi)*sin(psi) cos(theta)*cos(psi)];

end

%[appendix]{"version":"1.0"}
%---
