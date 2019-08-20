%% Cool fig for NATCOMM
% DKS 2019


%% config
config.path_dir = 'C:\Users\HE BEC\Dropbox\phd\projects\bell_witness\analysis\exp1_bell_max_viol\zxy_raw';
config.data_fname = 'exp1_1_raw-zxy_20190815';     % mat file

config.path_data = fullfile(config.path_dir,config.data_fname);


%% load data
S = load(config.path_data);

xyz_all = cellfun(@(x) zxy2xyz(x),S.zxy,'uni',0);


%% preproc
%%% config
prep_config.idx_collate = 1:10;    %50;

% cylindrical filter around halo centre to eliminate fringe hotspots
prep_config.det_cent=[0,3e-3,0];
prep_config.hotspot_rad=30e-3;

% XYZ ROI box 
prep_config.lim_x = 30e-3*[-1,1];
% prep_config.lim_y = [-1.2e-2,1.5e-2];
prep_config.lim_y = 30e-3*[-1,1];
prep_config.lim_z = [1.61,1.753];


%%% collate shots
xyz_collate = cat(1,xyz_all{prep_config.idx_collate});
n_shot_collate = numel(prep_config.idx_collate);


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


%% Mj distinguish
limz0 = prep_config.lim_z;      % original limz

%%% config
str_mJ = {'1','0'};
zsplit_mJ = mean(limz0);
cmap_mJ = {'blue','red'};


%
b_mJ0 = xyz_filt(:,3)>zsplit_mJ;
b_mJ1 = ~b_mJ0;
BB = {b_mJ1,b_mJ0};

xyz_mJ = cellfun(@(b) xyz_filt(b,:),BB,'uni',0);


%% VIS
%%% config
vis_config.FisUnits = 'centimeters';
vis_config.FigPos = [0,0,8,8];
vis_config.FigColor = 'k';

vis_config.Marker = 'o';
vis_config.SizeData = 3;
vis_config.MarkerFaceAlpha = 0.33;
vis_config.MarkerFaceColor = {[0,0,1], [1,0,0]};



% PLOT
H = figure('Name','dual_halo_3Dscatter');
H.Units = vis_config.FisUnits;
H.Position = vis_config.FigPos;

% H.Renderer = 'opengl';  % 'painters'
H.Renderer = 'painters';
H.InvertHardcopy = 'off';

set(H,'Color',vis_config.FigColor);

% P = scatter3(xyz_filt(:,1),xyz_filt(:,2),xyz_filt(:,3));

for ii=1:2
    txyz = xyz_mJ{ii};
    
    P = scatter3(txyz(:,1),txyz(:,2),txyz(:,3));
    
    set(P,'Marker',vis_config.Marker);
    set(P,'SizeData',vis_config.SizeData);
    set(P,'MarkerEdgeColor','none');
    set(P,'MarkerFaceColor',vis_config.MarkerFaceColor{ii});
    set(P,'MarkerFaceAlpha',vis_config.MarkerFaceAlpha);
    
    hold on;
end

axis equal;
axis off;