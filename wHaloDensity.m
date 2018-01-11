function [nn,az,el] = wHaloDensity(zxy,nAz,nEl,sig,lim)
% [NN, Z_AZ, Z_EL] = WHALODENSITY(ZXY, NAZ, NEL, SIG, LIM)
%
%

% define spherical zones
[az,el]=sphgrid(nAz,nEl);

% Cart to sph-polar
vs_sph=zxy2sphpol(zxy);
clearvars zxy;

% evaluate count density in sphere by the weighted counting algorithm
nCountsTot=size(vs_sph,1);
nn=zeros(size(az));
for ii=1:numel(az)
% parfor ii=1:numel(az)       % slower
    tzone=[az(ii),el(ii),1];
    tww=weightedCountSph(vs_sph,tzone,sig,lim);
    nn(ii)=tww;
end

% normalise weighted count density
wNormFactor=nCountsTot*normSphCount(nn,az,el);
nn=nn*wNormFactor;


end