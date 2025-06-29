%[text] # Equations of motion for orbit using Cowell
%[text] 軌道計算のための運動方程式
%[text] ## state variables
%[text] ${\\bf x}=\[x,y,z,v\_x,v\_y,v\_z\]^T$
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20210209  y.yoshimura, y.yoshimula@gmail.com
%[text] See also orbitConst.
function dxdt = eomOrbit(~, x_, const)
%[text] ## state variables and inputs
x = x_(1);
y = x_(2);
z = x_(3);
r = [x; y; z];

vx = x_(4);
vy = x_(5);
vz = x_(6);
v = [vx; vy; vz];

%[text] ## equations of motion
%[text] do not rewrite the following script
dxdt = [v
    -const.GEm ./ norm(r)^3 .* r ];

end


%[appendix]{"version":"1.0"}
%---
