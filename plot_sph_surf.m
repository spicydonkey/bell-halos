function ax=plot_sph_surf(azim,elev,ff)
% plots a 2D color-surface on sphere
% AX = PLOT_SPH_SURF(AZIM,ELEV,FF)
%
% AZIM: grid of azimuthal angles (rad); ordered and evenly distributed
% in [-pi,pi]; must be unique modulo 2*pi;
% ELEV: grid of elevation angles (meas from XY plane) (rad); ordered 
% and in range [-pi/2,pi/2]; must be unique;
% FF: grid of function values at (AZIM,ELEV) locations on sphere points
%
% ax: graphics object

% wrap angles to close surface
azim(end+1,:)=azim(1,:)+2*pi;
elev(end+1,:)=elev(1,:);
ff(end+1,:)=ff(1,:);

% sph-pol --> cartesian: can use surf!
r=1;    % mapped on a unit sphere
[xx,yy,zz]=sph2cart(azim,elev,r);

ax=surf(xx,yy,zz,ff);

set(ax,'FaceAlpha',0.7);

axis equal;
box off;
grid off;
axis off;
shading flat;
colormap magma;
