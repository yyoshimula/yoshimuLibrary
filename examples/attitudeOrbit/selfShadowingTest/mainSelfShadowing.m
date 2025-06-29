% ----------------------------------------------------------------------
%    self shadowing test using 3 patches (triangular patches)
%    20200811  y.yoshimura
%    Inputs:
%   Outputs:
%   related function files:
%   note:
%   cf: /lightcurves/staticLC/calcTest.mを改良
%   revisions;
%
%   (c) 2020 yasuhiro yoshimura
%----------------------------------------------------------------------

clc
clear
close all

global MAX_EQUIP EQUIP_NUM MAX_F
makeSat

sun = [-0.;0.0;1.0]; % nominal sun vector
sun = sun./norm(sun);
sun_dist = 1.496e+8 * 10^3;

face = SRP_LPS(face, sun', sun_dist);
face = calcShadow(face, vert, sun);

%% self shadowing
%
% 	// 平面の方程式を求める
% 	// 平面上の点P1(px1,py1,pz1), 点P2(px2,py2,pz2), 点P3(px3,py3,pz3)
% 	// 平面の方程式 a*x + b*y + c*z = d; ・・・(1)
% 	// また、P1P2 = u = (u1,u2,u3), P1P3 = v = (v1,v2,v3) とすると、
% 	// 平面の方程式で下記の式が成立。
% 	// |u2 v2| |u3 v3| |u1 v1|
% 	// |u3 v3| * (x-x0) + |u1 v1|*(y-y0) + |u2 v2|*(z-z0) = 0 ・・・（2）


%% show figures
figure(4460)
hold on
cb = colorbar;
colormap bone
% caxis([0 1e-7])

for m = 1:MAX_EQUIP
    for k = 1:EQUIP_NUM(m)
        h(m,k) = patch('Vertices', vert(m,k).coord, 'Faces', face(m,k).make);
        
        %         srpF = vecnorm(face(m,k).force,2,1); % norm of SRP force vector
        srp_ = face(m,k).flag .* -face(m,k).force(:,3);
        set(h(m,k), 'FaceColor', 'b');
        set(h(m,k), 'FaceLighting','phong','EdgeLighting','phong');
        set(h(m,k), 'FaceVertexCData', srp_, 'FaceColor', 'flat')
        axis equal;
        %     axis([-5,5, -5, 5, -5,5])
        grid on;
        xlabel('X'),ylabel('Y'),zlabel('Z')
        view(130,30);
        
    end
end
quiver3(0, 0, 0, 3*sun(1),3*sun(2),3*sun(3) ,'r')

%% variable attitude
atti = 0:50:50;

for i = 1:length(atti)
    DCM = ZYX2DCM(0.0, 0, deg2rad(atti(i)))'; % NOTE: vector rotation
    sun_var = DCM * sun;

    face = SRP_LPS(face, sun_var, sun_dist);
    face = calcShadow(face, vert, sun_var);

    quiver3(0, 0, 0, 3*sun_var(1),3*sun_var(2),3*sun_var(3) ,'r')

    for m = 1:MAX_EQUIP
        for k = 1:EQUIP_NUM(m)
            %             srpF = vecnorm(face(m,k).force,2,1); % norm of SRP force vector
            normedSRP = face(m,k).flag;% .* -face(m,k).force(:,3);
%             normedSRP = -face(m,k).force(3,:)';
            set(h(m,k), 'FaceVertexCData', normedSRP);%, 'FaceColor', 'flat');
        end
    end
    pause(0.5)
    
    movFrame(i) = getframe;

end

%% variable attitude with self-shadowing
% atti = 0:6:360;
% srp_fig = struct;
% for i = 1:length(atti)
%     DCM = ZYX2DCM(0.0, deg2rad(atti(i)), 0.0)'; % NOTE: vector rotation
%     sun_var = DCM * sun;
%     
%     face = SRP_LPS(face, sun_var, sun_dist);
%     face = calcShadow(face, vert, sun_var);
%     
%     for m = 1:MAX_EQUIP
%         for k = 1:EQUIP_NUM(m)
%             srp_fig(m,k).force(:,1) = face(m,k).flag .* -face(m,k).force(3,:)';
%         end
%     end
% end
% 
% for i = 1:length(atti)
%     DCM = ZYX2DCM(0.0, 0.0, deg2rad(atti(i)))'; % NOTE: vector rotation
%     sun_var = DCM * sun;
%     quiver3(0, 0, 0, 3*sun_var(1),3*sun_var(2),3*sun_var(3) ,'r')
%     
%     for m = 1:MAX_EQUIP
%         for k = 1:EQUIP_NUM(m)
%             set(h(m,k), 'FaceVertexCData', srp_fig(m,k).force, 'FaceColor', 'flat');
%         end
%     end
%     pause(0.4)
% end
