function [n,az,el]=plot_sphdist(kk,nBins)
% Simple diagnostic tool for plotting count distribution around the sphere
%
% [nn,az,el]=plot_sphdist(kk,n_az,n_el)
%
% kk:       nCounts x 3 vector array (zxy) of counts
% nBins:    1x2 - number of bins to uniformly divide azimuthal/elevation angles
%
% n:        count density in spherical zone
% az/el:    2D-grid of azimuthal/elevation angles to spherical zones
%

% parse inputs (defaults)
if ~exist('nBins','var')
    nBins=[100,50];
end

%%% density around sphere
% defaults
sig=[0.2,Inf];
lim=[3,Inf];

% count and normalise by gaussian convolution
[n,az,el]=haloZoneCount(kk,nBins(1),nBins(2),sig,lim,'gauss');


%%% Plot results
plot_sph_surf(az,el,n);

% annotate
cbar=colorbar;
cbar.Label.String='density [arb u]';

end