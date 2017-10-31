function [nn,z_az,z_el] = wHaloDensity(zxy,nAz,nEl,sig,lim)
% [NN, Z_AZ, Z_EL] = WHALODENSITY(ZXY, NAZ, NEL, SIG, LIM)
%
%

% define spherical zones
az=linspace(-pi,pi,nAz);
az=az(1:end-1);             % unique angles only
el=linspace(-pi/2,pi/2,nEl);

[z_az,z_el]=ndgrid(az,el);      % grid to preserve dimensional ordering in array indexing

% Cart to sph-polar
vs_sph=zxy2sphpol(zxy);
clearvars zxy;

% evaluate count density in sphere by the weighted counting algorithm
nCountsTot=size(vs_sph,1);
nn=zeros(size(z_az));
for ii=1:numel(z_az)
% parfor ii=1:numel(z_az)       % is slower
    tzone=[z_az(ii),z_el(ii),1];
    tww=weightedCountSph(vs_sph,tzone,sig,lim);
    nn(ii)=tww;
end

% normalise weighted count density
wNormFactor=nCountsTot*normSphCount(nn,z_az,z_el);
nn=nn*wNormFactor;


end