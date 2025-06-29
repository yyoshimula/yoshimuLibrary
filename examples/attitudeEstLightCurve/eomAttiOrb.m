%[text] # Equations of motion for attitude and orbit using Cowell
%[text] 姿勢・軌道計算のための運動方程式
%[text] ## state variables
%[text] ${\\bf x}=\[x,y,z,v\_x,v\_y,v\_z, {\\bf q}^T, {\\bf \\omega}^T\]^T$
%[text] `q`: quaternion where q(4) is the scalar part 
%[text] `w`: angular rate, rad/s
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20210209  y.yoshimura, y.yoshimula@gmail.com
%[text] See also orbitConst, mainATOM.
function dxdt = eomAttiOrb(~, x_, const, sat)
%[text] ## state variables and inputs
x = x_(1);
y = x_(2);
z = x_(3);
r = [x; y; z];

vx = x_(4);
vy = x_(5);
vz = x_(6);
v = [vx; vy; vz];

q1 = x_(7); %
q2 = x_(8);
q3 = x_(9);
q4 = x_(10);
q = [q1; q2; q3; q4];
wx = x_(11);
wy = x_(12);
wz = x_(13);
w = [wx; wy; wz];

tEx = zeros(3,1);
aEx = zeros(3,1);

MOI = sat.MOI;
  
kine = 0.5 .* qMult(4, 1, [w',0], q');
dwdt = MOI^(-1) * (cross(-w, MOI*w) + tEx);
%[text] ## equations of motion
%[text] do not rewrite the following script
dxdt = [v
    -const.GEm ./ norm(r)^3 .* r + aEx
    kine'
    dwdt   ];

end


%[appendix]{"version":"1.0"}
%---
