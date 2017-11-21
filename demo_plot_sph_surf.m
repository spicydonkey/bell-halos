% demo - plot_sph_surf

n_azim=20;
n_elev=20;
azim=linspace(0,2*pi,n_azim);
elev=linspace(-pi/2,pi/2,n_elev);
[AZIM,ELEV]=ndgrid(azim,elev);

nn=cos(ELEV)+sin(AZIM);

plot_sph_surf(AZIM,ELEV,nn);
