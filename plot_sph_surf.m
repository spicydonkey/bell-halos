function ax=plot_sph_surf(azim,elev,yy)
% plots a 2D color-surface on sphere
% plot_sph_surf(azim,elev,yy)
%
% azim: meshgrid of azimuthal angles (rad)
% elev: meshgrid of elevation angles (meas from XY plane) (rad)
% yy: meshgrid of double array of data at respective point on sphere
%
% ax: graphics object

% sph-pol --> cartesian: can use surf!
r=1;    % mapped on a unit sphere
[X,Y,Z]=sph2cart(azim,elev,r);

hs=surf(X,Y,Z,yy);

set(hs,'FaceAlpha',0.7);

axis equal;
box off;
grid off;
axis off;
shading flat;
colormap magma;