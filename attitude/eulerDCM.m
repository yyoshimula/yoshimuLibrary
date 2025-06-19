%[text] # Euler angle to directional cosine matrix (rotation matrix)
%[text] ## inputs
%[text] `phi`: rotation angle around the 1st axis (1st rotation)
%[text] `theta`: rotation angle around 2nd axis (2nd rotation)
%[text] `psi`: rotation angle around 3rd axis (3rd rotation)
%[text] ## output
%[text] R: directional cosine matrix (rotation matrix), 3x3
%[text] ## note
%[text] if rotation matrix of ZYX (3-2-1) Euler angle is required:
%[text] ```matlabCodeExample
%[text] phi = deg2rad(30);
%[text] theta = deg2rad(20);
%[text] psi = deg2rad(10);
%[text] 
%[text] R = eulerDCM(3, 2, 1, phi, theta, psi);
%[text] ```
%[text] ## references 
%[text] Markley, F. L., and Crassidis, J. L., Fundamentals of Spacecraft Attitude Determination and Control, New York, NY: Springer, 2014. 
%[text] ## revisions
%[text] 20150101  y.yoshimura, y.yoshimula@gmail.com
%[text] See also q2DCM.
function R = eulerDCM(axis1, axis2, axis3, phi, theta, psi)
% arguments
%     axis1 (1,1) {mustBeMember(axis1, [1, 2, 3])}
%     axis2 (1,1) {mustBeMember(axis2, [1, 2, 3])}
%     axis3 (1,1) {mustBeMember(axis3, [1, 2, 3])}    
%     phi (1,1) 
%     theta (1,1)
%     psi (1,1) 
% end

R = dcm1axis(axis3, psi) * dcm1axis(axis2, theta) * dcm1axis(axis1, phi);

end

%[appendix]{"version":"1.0"}
%---
