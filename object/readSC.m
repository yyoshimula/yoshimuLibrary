%[text] # Reading spacecraft shape data 
%[text] 衛星形状や表面特性の読み込み
%[text] `satName` : satellite object file name, .obj file
%[text] `sat`: 構造体satに各種データを格納
%[text] `sat.vertices`: vertex position, (x, y, z), m, Nx3 matrix, 面を構成する点の座標
%[text] `sat.normal`: normal vector, Nx3 matrix, 法線方向単位ベクトル, 外向きが正
%[text] `sat.faces`: face indices, Nx3 matrix, 面を構成する座標のindex
%[text] `sat.area`: face area, ${\\rm m^2}$, Nx1 vector, 面の面積
%[text] `sat.pos`: center of face, m, Nx3 matrix, 面の中心（平均）座標
%[text] `sat.Ca`: coefficients for absorption, Nx1 vector, 吸収率
%[text] `sat.Cd`: coefficients for diffusion, Nx1 vector, 拡散反射率
%[text] `sat.Cs`: coefficients for specular reflection, Nx1 vector, 鏡面反射率
%[text] `sat.sunlitFlag`: self shadow flag, (default)1: not shadowoed, 0: shadowed
%[text] ## note
%[text] ## references
%[text] NA
%[text] ## revisions
%[text] 20240729 major update, readObjで読み込むように変更 y.yoshimura
%[text] 2021020209  y.yoshimura
%[text] See also showSC.
function sat = readSC(satName)

sat = readObj(satName);

[sat.area, sat.pos] = calcAreaObj(sat);  % face areas and position

[sat.uu, sat.uv, sat.qlb] = calcLocalFrame(sat); % local frames and their quaternions

N = size(sat.area,1);
%[text] ### initialize
% reflectivity
sat.F0 = ones(N,1) .* 0.5; % default value
sat.kappa = zeros(N,1);
sat.Ca = zeros(N,1);
sat.Cd = ones(N,1) .* 0.5; % default value
sat.Cs = ones(N,1) .* 0.5; % default value
sat.nu = ones(N,1) .* 800; % default value
sat.nv = ones(N,1) .* 800; % default value
sat.fObs = zeros(N,1);

% physical (default values)
sat.MOI = diag([10, 15, 20]);
sat.m = 100;

sat.sunlitFlag = ones(N,1);

sat.forces = zeros(N,3);
sat.torque = zeros(N,3);

disp('number of vertices:')
disp(size(sat.vertices, 1)) % nVertices
disp('number of faces:')
disp(size(sat.faces, 1)) % nFaces

end

%[appendix]{"version":"1.0"}
%---
