function [r0,dr,az,el]=summary_disthalo_ellip(K,nBins,dtheta,verbose)
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
if ~exist('verbose','var')
    verbose=1;
end

% preproc
k_all=collate_shots(K);     % we don't need corr info
[az,el]=sphgrid(nBins(1),nBins(2));
nspecies=size(K,2);

% peak radii around halo
r0=cell(1,nspecies);
dr=cell(1,nspecies);
parfor ii=1:nspecies
    [r0{ii},dr{ii}]=arrayfun(@(a,e) localRadDist(double(k_all{ii}),...
        a,e,dtheta),az,el);
end
    
% [r0,dr]=cellfun(@(c)...
%     arrayfun(@(a,e) localRadDist(double(c),a,e,dtheta,0),az,el),...
%     k_all,'UniformOutput',false);


%% plots
if verbose>0
    for ii=1:nspecies
        figure;
        % radial distribution
        subplot(1,2,1);
        plotFlatMapWrappedRad(az,el,r0{ii},'eckert4');
        cbar=colorbar('southoutside');
        cbar.Label.String='Radius';
    
        % width distribution
        subplot(1,2,2);
        plotFlatMapWrappedRad(az,el,dr{ii},'eckert4');
        cbar=colorbar('southoutside');
        cbar.Label.String='rms width';
    end
end

end