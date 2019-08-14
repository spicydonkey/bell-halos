%% Create 2D images of Bell data (publication fig.1)
% DKS
% 20180822
%
% Prereq: prepare Nx3 "xyz"-vector array of data


%% CONFIG
% path do data directory
% path_dir='C:\Users\HE BEC\Documents\lab\bell_momentumspin\bell_epr_2018\proc\fig1_halo_img';
config.path_dir='C:\Users\David\Dropbox\PhD\projects\bell_witness\paper\fig1\fig1_halo_img';

% data filename
% config.data_fname = 'theta_0.mat';         % no rotation (SRC)
config.data_fname = 'theta_pi2.mat';     % pi/2 rotation (ROT): too many shots


%% load data
config.path_data=fullfile(config.path_dir,config.data_fname);   

S = load(config.path_data);


%% preprocess
% TODO: collate shots and save #shots for normalisation

%%% config
% cylindrical filter around halo centre to eliminate fringe hotspots
prep_config.det_cent=[0,3e-3,0];
prep_config.hotspot_rad=30e-3;

% XYZ ROI box 
prep_config.lim_x = 30e-3*[-1,1];
prep_config.lim_y = [-1.2e-2,1.5e-2];
prep_config.lim_z = [1.61,1.753];


%%% filter data
xyz_filt=S.xyz;

% cylindrical filter
xyz_filt = cylindercull(xyz_filt,prep_config.det_cent,[prep_config.hotspot_rad,Inf],3);

% box filter
xyz_filt = boxcull(xyz_filt,{prep_config.lim_x,prep_config.lim_y,prep_config.lim_z});


%%% VIS
H_preproc=figure('Name','scatter_preproc');
scatter3(xyz_filt(:,1),xyz_filt(:,2),xyz_filt(:,3),'b.');
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');


%% 2D counts histogram
%%% config
hist_config.voxel = 8e-4*[1,1,1];     % voxel size along x,y,z-axis
hist_config.lim = {prep_config.lim_x,prep_config.lim_y,prep_config.lim_z};   % histogram limits


hist_config.ed = cellfun(@(x,dx) x(1):dx:x(2),hist_config.lim,...
    num2cell(hist_config.voxel),'uni',0);   % bin edges
hist_config.Vvoxel = prod(hist_config.voxel);   % voxel volume in m^3 


%%% histogram
% 3D hist
[count_3d,edges,mid] = histcn(xyz_filt,hist_config.ed{1},hist_config.ed{2},hist_config.ed{3});     
countden_3d = count_3d/(hist_config.Vvoxel*1e9);      % count number density (in mm^-3)

% 2D projection: XZ-plane
countden_xz = squeeze(mean(countden_3d,2));


%%% vis
H_2d = figure('Name','countdensity');
Im = imagesc(mid{1}, mid{3}, log(countden_xz)');

cbar = colorbar;
cbar.Label.String = 'atom density (mm^{-3})';
cbar.TickLabels = cellstr(num2str(cbar.Ticks(:), '10^%d'));
% XTickLabels = cellstr(num2str(round(log10(XTick(:))), '10^%d'));
colormap('blue');

set(gca,'YDir','Normal');
axis equal;
axis tight;

