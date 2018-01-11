function [h,nn,az,el]=summary_disthalo_ndist(K,nBins,plotmode)
% Summarise distinguishable (2-species) halos with density distribution
%
% [h,nn,az,el]=summary_disthalo_ndist(K,OPTS)
%

% check args
if ~exist('nBins','var')
    warning('nBins is undefined. Setting to default value: [100,50].');
    nBins=[100,50];
end
if ~exist('plotmode','var')
    warning('plotmode is undefined. Setting to default value: flat.');
    plotmode='flat';
end

% configs
color_mf={'b','r'};

%%% radial dist
h(1)=figure();
for ii=1:2
    hold on;
    [~,~,tbar]=plot_rdist(vertcat(K{:,ii}));
    tbar.FaceColor=color_mf{ii};
end
drawnow

%%% sph density dist
h(2)=figure();
[az,el]=sphgrid(nBins(1),nBins(2));
nn=cell(1,2);       % number density (arb) around halo
for ii=1:2
    subplot(1,2,ii);
    [nn{ii}]=plot_sphdist(vertcat(K{:,ii}),nBins,plotmode);
end
drawnow

end