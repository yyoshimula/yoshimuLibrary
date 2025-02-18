%[text] # Visualizing spacecraft
%[text] 読み込んだ衛星形状を表示
%[text] `sat`: 各種データが格納された構造体
%[text] `sat.vertices`: vertex position, (x, y, z), m, Nx3 matrix, 面を構成する点の座標
%[text] `sat.faces`: face indices, Nx3 matrix, 面を構成する座標のindex
%[text] `sat.area`: face area, ${\\rm m^2}$, Nx1 vector, 面の面積
%[text] `sat.pos`: center of face, m, Nx3 matrix, 面の中心（平均）座標
%[text] `sat.Ca`: coefficients for absorption, Nx1 vector, 吸収率
%[text] `sat.Cd`: coefficients for diffusion, Nx1 vector, 拡散反射率
%[text] `sat.Cs`: coefficients for specular reflection, Nx1 vector, 鏡面反射率
%[text] `sat.n`: normal vector, Nx3 matrix, 法線方向単位ベクトル, 外向きが正
%[text] `sat.shadowFlag`: self shadow flag, (default)1: not shadowoed, 0: shadowed
%[text] `options`
%[text] `    Normal`: draw normal vectors of facets, `on` or `off`
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 2021020127  y.yoshimura
%[text] See also readSC.
function showSC(sat, options)
arguments
    sat
    options.Normal char {mustBeMember(options.Normal,{'on', 'off'})} = 'off' 
end

figure
patch('Faces',sat.faces,'Vertices',sat.vertices, 'facecolor', [0.7 0.71 0.71], 'facealpha', 0.9);
if strcmp(options.Normal,'on') % normal vector arrows
    hold on
    quiver3(sat.pos(:,1), sat.pos(:,2), sat.pos(:,3), sat.normal(:,1), sat.normal(:,2), sat.normal(:,3));
    quiver3(sat.pos(:,1), sat.pos(:,2), sat.pos(:,3), sat.uu(:,1), sat.uu(:,2), sat.uu(:,3)); % x-axis of local frame
end


axis equal
xlabel('x'), ylabel('y'), zlabel('z')
view([28.0 42.0])
grid on

end

%[appendix]{"version":"1.0"}
%---
