%% halo density images for experiment schematic (Fig.1)
% DKS 2019

%% config
config.path_dir = 'C:\Users\HE BEC\Dropbox\phd\projects\bell_witness\analysis\exp1_bell_max_viol\zxy_raw';
config.data_fname = 'exp1_1_raw-zxy_20190815';     % mat file

config.path_data = fullfile(config.path_dir,config.data_fname);


%% load data
S = load(config.path_data);

xyz_all = cellfun(@(x) zxy2xyz(x),S.zxy,'uni',0);
 

%% preprocess
%%% config
% shots
idx_shot_collate = 1:1000;     % shot indices to collate

% cylindrical filter around halo centre to eliminate fringe hotspots
prep_config.det_cent=[0,3e-3,0];
prep_config.hotspot_rad=30e-3;

% XYZ ROI box 
prep_config.lim_x = 30e-3*[-1,1];
prep_config.lim_y = [-1.2e-2,1.5e-2];
% prep_config.lim_y = 10e-3*[-1,1];
prep_config.lim_z = [1.61,1.753];


%%% collate shots
xyz_collate = cat(1,xyz_all{idx_shot_collate});
n_shot_collate = numel(idx_shot_collate);


%%% filter data
xyz_filt=xyz_collate;

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
% histogram
hist_config.voxel = 8e-4*[1,1,1];     % voxel size along x,y,z-axis
hist_config.lim = {prep_config.lim_x,prep_config.lim_y,prep_config.lim_z};   % histogram limits


%%% histogram
hist_config.ed = cellfun(@(x,dx) x(1):dx:x(2),hist_config.lim,...
    num2cell(hist_config.voxel),'uni',0);   % bin edges
hist_config.Vvoxel = prod(hist_config.voxel);   % voxel volume in m^3 


% 3D hist
[count_3d,edges,mid] = histcn(xyz_filt,hist_config.ed{1},hist_config.ed{2},hist_config.ed{3});     
countden_3d = count_3d/((hist_config.Vvoxel*1e9) * n_shot_collate);      % avg single-shot count number density (in mm^-3)

% 2D projection: XZ-plane
countden_xz = squeeze(mean(countden_3d,2));


%% VIS all mJ
%%% config
% width = 3.3;
% height = 5.6/2;
% fontsize = 9;

% scale = h.Position(3)/width;        % to scale width/height by

%%%
H = figure('Name','countdensity');
% H.Units = 'centimeters';
% H.Position = [0,0,width,height];

Im = imagesc(mid{1}, mid{3}, log(countden_xz)');

cbar = colorbar;
cbar.Label.String = 'atom density (mm^{-3})';
cbar.TickLabels = cellstr(num2str(cbar.Ticks(:), '10^%d'));
colormap('red');

set(gca,'YDir','Normal');
axis equal;
axis tight;



%% vis
limz0 = prep_config.lim_z;      % original limz

%%% config
str_mJ = {'1','0'};
zsplit_mJ = mean(limz0);
cmap_mJ = {'blue','red'};

fontsize = 17.5;

%%% imgs
limz_mJ = {[limz0(1),zsplit_mJ],[zsplit_mJ,limz0(2)]};

H_2d=[];
Im=[];
for ii=1:2
    H_2d(ii) = figure('Name',['countdensity_',str_mJ{ii}]);
%     H.Units = 'centimeters';
%     H.Position = [0,0,width,height];
    
    Im(ii) = imagesc(mid{1}, mid{3}, log(countden_xz)');
    
    colormap(cmap_mJ{ii});
    
    cbar = colorbar;
    cbar.Label.String = 'avg. count density (mm^{-3})';
    cbar.Ticks = cbar.Ticks(1:4:end);
    cbar.TickLabels = cellstr(num2str(cbar.Ticks(:), '10^{%d}'));
%     if ii == 1
%         cbar.Direction='reverse';
%     end
    
    %
    h = gcf;
    ax = gca;
    
    set(gca,'YDir','Normal');
    axis equal;
    
    xlim(prep_config.lim_x);
    ylim(limz_mJ{ii});

    axis off;
    
    set(findall(h,'-property','FontSize'),'FontSize',fontsize);
end
