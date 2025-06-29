% ----------------------------------------------------------------------
%    self shadowing test using 3 patches (triangular patches)
%    20200811  y.yoshimura
%    Inputs:
%   Outputs:
%   related function files:
%   note:
%   cf: /lightcurves/staticLC/calcTest.m‚ð‰ü—Ç
%   revisions;
%
%   (c) 2020 yasuhiro yoshimura
%----------------------------------------------------------------------

clc
clear
close all

global MAX_EQUIP EQUIP_NUM MAX_F
makeSat

sun = [0;0.0;1.0]; % nominal sun vector
sun = sun./norm(sun);
sun_dist = 1.496e+8 * 10^3;

face = SRP_LPS(face, sun, sun_dist);
face = calcShadow(face, vert, sun);

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
        srp_ = face(m,k).flag .* -face(m,k).force(3,:)';
        set(h(m,k), 'FaceColor', 'b');
        set(h(m,k), 'FaceLighting','phong','EdgeLighting','phong');
        set(h(m,k), 'FaceVertexCData', srp_, 'FaceColor', 'flat')
        axis vis3d equal;
        %     axis([-5,5, -5, 5, -5,5])
        grid on;
        xlabel('X'),ylabel('Y'),zlabel('Z')
        view(130,30);
        
    end
end
quiver3(0, 0, 0, 3*sun(1),3*sun(2),3*sun(3) ,'r')

%% csv output
outData = [];
outData2 = [];
for i = 1:MAX_EQUIP
    for j = 1:EQUIP_NUM(i)        
        outData = [outData
            face(i,j).area', face(i,j).normal', ...
            face(i,j).Ca', face(i,j).Cd', face(i,j).Cs', face(i,j).Ct',...
            face(i,j).make];
        
         outData2 = [outData2
            vert(i,j).coord];
           
    end
end

writematrix(outData,'face.txt','Delimiter', ' ');
writematrix(outData2,'vert.txt','Delimiter', ' ');

% % csvwrite('face.csv', outData);
% csvwrite('vert.csv', outData2);