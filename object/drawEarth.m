%[text] # visualizing Earth
%[text] ## note
%[text] NA
%[text] ## references 
%[text] ## revisions
%[text] 20220612  y.yoshimura, y.yoshimula@gmail.com
%[text] See also orbitConst.
function drawEarth(GMST, alpha, const)

backColor = 'w';
nPanels = 180;

% Earth texture image
image_file = 'earth.jpg';

% Mean spherical earth

fig = gcf;
fig.Color = backColor;
hold on;

[x, y, z] = ellipsoid(0, 0, 0, const.RE, const.RE, const.RE, nPanels);
globe = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);

hgx = hgtransform;
tmp = makehgtform('zrotate',GMST);
set(hgx,'Matrix', tmp);
set(globe,'Parent', hgx);

%[text] ## texture mapping
cdata = imread(image_file);

set(globe,'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');
drawnow

end

%[appendix]{"version":"1.0"}
%---
