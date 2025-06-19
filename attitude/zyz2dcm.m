%[text] # ZYZ (3-2-3) Euler angle to directional cosine matrix (rotation matrix)
%[text] ## inputs
%[text] `phi`: rotation angle around Z-axis (1st rotation)
%[text] `theta`: rotation angle around Y-axis (2nd rotation)
%[text] `psi`: rotation angle around Z-axis (3rd rotation)
%[text] ## output
%[text] R: directional cosine matrix (rotation matrix), 3x3
%[text] ## note
%[text] NA
%[text] ## references 
%[text] Markley, F. L., and Crassidis, J. L., Fundamentals of Spacecraft Attitude Determination and Control, New York, NY: Springer, 2014. 
%[text] ## revisions
%[text] 20150101  y.yoshimura, y.yoshimula@gmail.com
%[text] See also q2DCM.
function R = zyz2dcm(phi, theta, psi)
% arguments
%     phi (1,1) 
%     theta (1,1)
%     psi (1,1) 
% end

R = dcm1axis(3, psi) * dcm1axis(2, theta) * dcm1axis(3, phi);

end

%[appendix]{"version":"1.0"}
%---
