% Anomalous Rabi-flopping
%   Investigate non-flopping population oscillations in some regions of the
%   halo under a single-beam Raman process.
%
% DK Shin
%


% Load data
load('theta_rabi_20180131.mat');


% get sph-grid
n_sphgrid=size(az);
n_az=n_sphgrid(1);
n_el=n_sphgrid(2);


% select parts of halo sph-zones
el_in=find(abs(el(1,:))<0.8);         % polar coords away from poles

idx_az=1:13:n_az;
idx_el=el_in(1):5:el_in(end);


% Plot oscillations
% TODO - do for mf=0 data after background subtraction
this_mf_idx=2;      % 1 for mf=0; 2 for mf=1 

data_color=distinguishable_colors(numel(idx_az));          % thetha
data_linestyle={'-','--',':','-.'};      % phi

figure('Name','rabi_cycle');
for ii=1:numel(idx_az)
    for jj=1:numel(idx_el)
        hold on;
        plot(1e6*traman,squeeze(P_mf{this_mf_idx}(idx_az(ii),idx_el(jj),:)),...
            'Color',data_color(ii,:),...
            'LineStyle',data_linestyle{jj},...
            'LineWidth',1.5);
        
    end
end
xlabel('$\tau$ ($\mu$s)');
ylabel('$P$');
box on;