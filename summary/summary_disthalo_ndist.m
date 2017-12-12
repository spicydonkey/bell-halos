function [h]=summary_disthalo_ndist(K)
% Summarise distinguishable (2-species) halos with density distribution
%
% [h]=summary_disthalo_ndist(K)
%

%%% configs
color_mf={'b','r'};

% radial dist
h(1)=figure();
for ii=1:2
    hold on;
    [~,~,tbar]=plot_rdist(vertcat(K{:,ii}));
    tbar.FaceColor=color_mf{ii};
end
drawnow

% sph dist
h(2)=figure();
for ii=1:2
    subplot(1,2,ii);
    plot_sphdist(vertcat(K{:,ii}));
end
drawnow

end