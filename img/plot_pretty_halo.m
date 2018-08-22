function [n_v,v_ed,xyz_filt] = plot_pretty_halo(xyz,dbg)
% making pretty halo images
%   

if ~exist('dbg','var')
    dbg=false;
end

xyz_raw=xyz;


%% filter raw XYZ
% filter detector hot-spots
%   cylinrical filter
det_cent=[0,3e-3,0];        
hotspot_rad=30e-3;

[xyz_filt] = cylindercull(xyz_raw,det_cent,[hotspot_rad,Inf],3);


% filter to Z-window
lim_z=[1.61,1.753];
[xyz_filt] = boxcull(xyz_filt,{[],[],lim_z});


% filter to Y-slice around center of halo
lim_y=[-1.2e-2,1.5e-2];
[xyz_filt] = boxcull(xyz_filt,{[],lim_y,[]});

% vis: scatter plot
hfig=figure('Name','pretty_halo');
scatter3(xyz_filt(:,1),xyz_filt(:,2),xyz_filt(:,3),'b.');
axis equal;


%% to 2D density
n_z=4000;
range_xyz=range(xyz_filt);

n_bins=[round(n_z*range_xyz(1)/range_xyz(3)),1,n_z];
[n_v,v_ed]=histVectors(xyz_filt,n_bins);

n_v=squeeze(n_v)';      % squeeze to 2D (XZ)

% filter
gaussfilt_sig=10;
n_filt=imgaussfilt(n_v,gaussfilt_sig);
n_disp=log(n_filt);     % to display in arb unit
n_lim=[-6 -1];         	% Caxis limits hardset

% normalise displayable val to [0,1]
n_disp=(n_disp-n_lim(1))/diff(n_lim);

% vis: density image
hfig_2d=figure('Name','density');
imagesc(v_ed{1},v_ed{3},n_disp);
set(gca,'YDir','normal')
axis equal;
cbar=colorbar('southoutside','Ticks',0:0.5:1.0);
cbar.Label.String='Scaled density (arb. unit)';
cbar.FontSize=20;
caxis([0,1]);
colormap('inferno');


%% DEBUG
if dbg
    h_dbg=figure('Name','debug');
    
    % plot original point cloud
    scatter3(xyz_raw(:,1),xyz_raw(:,2),xyz_raw(:,3),'b.');
    hold on;
    
    % capt'd counts
    scatter3(xyz_filt(:,1),xyz_filt(:,2),xyz_filt(:,3),'ro');
end
end