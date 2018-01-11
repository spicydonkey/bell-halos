function [n,az,el]=plot_sphdist(kk,nBins,plotmode)
% Simple diagnostic tool for plotting count distribution around the sphere
%
% [nn,az,el]=plot_sphdist(kk,n_az,n_el)
%
% kk:       nCounts x 3 vector array (zxy) of counts
% nBins:    1x2 - number of bins to uniformly divide azimuthal/elevation angles
% plotmode: string for plotting type flat/3d
%
% n:        count density in spherical zone
% az/el:    2D-grid of azimuthal/elevation angles to spherical zones
%

% parse inputs (defaults)
if ~exist('nBins','var')
    warning('nBins is undefined. Setting to default value: [100,50].');
    nBins=[100,50];
end
if ~exist('plotmode','var')
    warning('plotmode is undefined. Setting to default value: flat.');
    plotmode='flat';
end
    
%%% density around sphere
% defaults
count_mode='gauss';
sig=[0.1,Inf];
lim=[3,Inf];

[az,el]=sphgrid(nBins(1),nBins(2));     % create grid
n=haloZoneCount(kk,az,el,sig,lim,count_mode);   % count around sphere

%%% Plot results
switch plotmode
    case 'flat'
        plotFlatMapWrappedRad(az,el,n,'eckert4');
    case '3d'
        plot_sph_surf(az,el,n);
    otherwise
        error('Unknown plotmode (flat/3d).');
end

% annotate
hcb=colorbar('Southoutside');
hcb.TickLabelInterpreter='latex';
hcb.Label.Interpreter='latex';
hcb.Label.String='density [arb u]';

% configs summary
%   nBins; sig; lim
titlestr=sprintf('[%d,%d]; [%0.3g,%0.3g]; [%0.3g,%0.3g]',nBins,sig,lim);
title(titlestr);

end