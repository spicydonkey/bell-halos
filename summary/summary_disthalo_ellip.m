function [r0,dr,az,el]=summary_disthalo_ellip(K,nBins,dtheta)
%Summarise ellipticity of halo
%   * Peak-radial distribution around halo
%
%   K: NSHOTxNSPEC cell-array of halos
%
%   r0: 1xNSPEC cell - peak radius
%   dr: 1xNSPEC cell - rms width of radial profile
%   az: azim grid
%   el: elev grid
%

% check args and config
if ~exist('nBins','var')
    warning('nBins is undefined. Setting to default value: [30,15].');
    nBins=[30,15];
end
if ~exist('dtheta','var')
    warning('dtheta is undefined. Setting to default value: 0.3.');
    dtheta=0.3;
end

% preproc
k_all=collate_shots(K);     % we don't need corr info
[az,el]=sphgrid(nBins(1),nBins(2));

% peak radii around halo
[r0,dr]=cellfun(@(c)...
    arrayfun(@(a,e) localRadDist(double(c),a,e,dtheta,0),az,el),...
    k_all,'UniformOutput',false);

% plots
figure;
for ii=1:numel(r0)
    subplot(1,numel(r0),ii);
    plotFlatMapWrappedRad(az,el,r0{ii},'eckert4');
    cbar=colorbar('southoutside');
    cbar.Label.String='Mean Radius';
end

end