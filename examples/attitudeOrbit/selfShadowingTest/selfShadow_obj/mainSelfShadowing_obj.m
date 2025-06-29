% ----------------------------------------------------------------------
%    self shadowing test using OBJ file reader
%    20200811  y.yoshimura
%    Inputs:
%   Outputs:
%   related function files:
%   note: OBJ file reader is required:
%         https://jp.mathworks.com/matlabcentral/fileexchange/27982-wavefront-obj-toolbox?s_tid=FX_rc1_behav
%   revisions;
% 
%   ----Structure Format-----
%   sat.vertices: vertex position, (x, y, z), m, Nx3 matrix, 面を構成する点の座標
%   sat.faces: face indices, Nx3 matrix, 面を構成する座標のindex
%   sat.area: face area, m^2, Nx1 vector, 面の面積
%   sat.pos: center of face, m, Nx3 matrix, 面の中心（平均）座標
%   sat.Ca: coefficients for absorption, Nx1 vector, 吸収率
%   sat.Cd: coefficients for diffusion, Nx1 vector, 拡散反射率
%   sat.Cs: coefficients for specular reflection, Nx1 vector, 鏡面反射率
%   sat.n: normal vector, Nx3 matrix, 法線方向単位ベクトル, 外向きが正
%   sat.shadowFlag: self shadow flag, (default)1: not shadowoed, 0: shadowed
%   (c) 2020 yasuhiro yoshimura
%----------------------------------------------------------------------

clc
clear
close all

%% load object
% いったん変数objに読み込む
obj = read_wobj('oneweb.obj');

nCompo = 3; % the number of components, modellingによって変更する
sat.vertices = obj.vertices; % all vertices, Nx3 matrix
tmp = [];
for i = 1:nCompo
    tmp = [tmp
        obj.objects(4*i).data.vertices];
end
sat.faces = tmp; % face indices
sat = calcArea_obj(sat);  % face areas
sat = calcNormal_obj(sat); % face normal vectors
sat.Ca = ones(length(sat.vertices),1) .* obj.material(3).data(1); % objファイルの定義はcolorだが，とりあえず
sat.Cd = ones(length(sat.vertices),1) .* obj.material(4).data(1);
sat.Cs = ones(length(sat.vertices),1) .* obj.material(5).data(1);

%% show object shape (for debug)
figure
patch('Faces',sat.faces,'Vertices',sat.vertices, 'facecolor', [0 0 1]);
hold on
quiver3(sat.pos(:,1), sat.pos(:,2), sat.pos(:,3), ...
    3*sat.normal(:,1), 3*sat.normal(:,2), 3*sat.normal(:,3),'r');
axis equal

%% parameters
sun = [0.8; 0.4; 0.0]; % nominal sun vector
sun = sun./norm(sun);
sun_dist = 1.496e+8 * 10^3; % m

sat = SRP_LPS_obj(sat, sun', sun_dist); % calc. SRP
sat = calcShadow_obj(sat, sun); % calc self shadowing

%% show figures
figure(4460)
hold on
cb = colorbar;
colormap bone
% caxis([0 1e-7])

satFig = patch('Faces',sat.faces,'Vertices',sat.vertices);
set(satFig, 'FaceColor', 'b');
set(satFig, 'FaceLighting','phong','EdgeLighting','phong');
set(satFig, 'FaceVertexCData', sat.shadowFlag, 'FaceColor', 'flat')
axis equal
grid on
xlabel('X'),ylabel('Y'),zlabel('Z')
view(130,30);
quiver3(0, 0, 0, 3*sun(1),3*sun(2),3*sun(3) ,'r') % show sun direction

%% calculation for variable attitude
atti = 0:-10:-80; % attitude angle , deg

for i = 1:length(atti)
    DCM = ZYX2DCM(deg2rad(atti(i)), 0.0, 0.0)'; % NOTE: vector rotation
    sun_var = DCM * sun;
    
    sat = SRP_LPS_obj(sat, sun_var, sun_dist); % SRP calculation
    sat = calcShadow_obj(sat, sun_var); % shadow flag calculation
    
    quiver3(0, 0, 0, 3*sun_var(1),3*sun_var(2),3*sun_var(3) ,'r')
    
    set(satFig, 'FaceVertexCData', sat.shadowFlag);
    
    pause(0.5)
    
%     movFrame(i) = getframe;
    
end

% movie(movFrame);