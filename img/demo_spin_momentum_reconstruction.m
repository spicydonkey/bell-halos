%% Halo (spin x momentum) space reconstruction
% prototyped for press-release
% DKS 2019



%% data for z-axis measurement (theta = 0)
% Analysis of global rotation on correlated atom-pair 
%
% 2018.05.04: raw basis-dependent spin correlation - low mode occupancy
%
% DKS
%

%%% EXP1 pi/2 rotation
% config_name='C:\Users\HE BEC\Documents\MATLAB\bell-halos\analysis\exp1_bell_max_viol\src\config_1.m';
% config_name='C:\Users\HE BEC\Documents\MATLAB\bell-halos\analysis\exp1_bell_max_viol\src\config_3.m';
% config_name='C:\Users\HE BEC\Documents\MATLAB\bell-halos\analysis\exp1_bell_max_viol\src\config_5.m';
% config_name='C:\Users\HE BEC\Documents\MATLAB\bell-halos\analysis\exp1_bell_max_viol\src\config_7.m';
% config_name='C:\Users\HE BEC\Documents\MATLAB\bell-halos\analysis\exp1_bell_max_viol\src\config_8.m';

%%% EXP5 Bell Y-rotation
config_name='C:\Users\HE BEC\Documents\MATLAB\bell-halos\analysis\exp5_bell_yrot\src\config_1.m';
% config_name='C:\Users\HE BEC\Documents\MATLAB\bell-halos\analysis\exp5_bell_yrot\src\config_2.m';
% config_name='C:\Users\HE BEC\Documents\MATLAB\bell-halos\analysis\exp5_bell_yrot\src\config_3.m';

%%% pi/4 Y-rotation (higher mocc)
% config_name='C:\Users\HE BEC\Documents\MATLAB\bell-halos\analysis\ideal_global_rotation_v4\src\config_1.m';


%%% EXPX: X-rotation <JyJy>
% config_name='C:\Users\HE BEC\Documents\MATLAB\bell-halos\analysis\exp5_bell_yrot\src\config_x1.m';
% config_name='C:\Users\HE BEC\Documents\MATLAB\bell-halos\analysis\exp5_bell_yrot\src\config_x2.m';
% config_name='C:\Users\HE BEC\Documents\MATLAB\bell-halos\analysis\exp5_bell_yrot\src\config_x3.m';


% config_name='C:\Users\David\Documents\MATLAB\bell-halos\analysis\exp5_bell_yrot\src\config_2.m';
% config_name='C:\Users\David\Documents\MATLAB\bell-halos\analysis\exp5_bell_yrot\src\config_1.m';

%% load config
run(config_name);


%% Get experimental params
if configs.flag.param_scan
    % load wfmgen log
    [params,id_in_param,param_id,Ipar]=wfmgen_log_parser(configs.path.paramlog);
    nparam=size(params,1);      % number of unique param-set
    
    % get searched param
    par_T=params;       % scanned pulse duration [s]
else    
    % defaults autoscan-related vars
    nparam=1;       % 1- since we don't search any params
    id_in_param={configs.load.id};    % all IDS to load
    
%     par_T=NaN;      % default param to NaN
    par_T=configs.misc.param;       % const param hardcoded in configs
end

%% load txy
[txy,fout]=load_txy(configs.load.path,configs.load.id,configs.load.window,configs.load.mincount,configs.load.maxcount,[],1,2,0);

zxy=txy2zxy(txy,vz);
shot_id=fout.id_ok;

% build boolean array selector for scanned param-set
b_paramset=cellfun(@(idx) ismember(shot_id,idx),...
    id_in_param,'UniformOutput',false);
b_paramset=horzcat(b_paramset{:});


% DEBUG
h_zxy_raw=figure;
plot_zxy(zxy,1e5);
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
view([0,0]);


%% distinguish mF and capture halo
n_shot=size(zxy,1);
n_mf=numel(configs.mf);

% preallocate
p_bec=cell(1,n_mf);     % marker BEC center positions
p_bec0=cell(1,n_mf);    % marker BEC in halo centered ZXY coord sys
p_halo=cell(1,n_mf);   % approx halo center (mid-point of marker BECs)
for ii=1:n_mf
    p_bec{ii}=NaN(n_shot,3,2);
    p_halo{ii}=NaN(n_shot,3);
end
zxy0=cell(n_shot,n_mf); % halo centerised zxy


for ii=1:n_shot
    tzxy=zxy{ii};
    for jj=1:n_mf
        tzxy_mf=boxcull(tzxy,configs.mf(jj).window);
        tp_bec0=configs.mf(jj).p_bec0;
        tr_bec0=configs.mf(jj).r_bec0;
        
        for kk=1:size(tp_bec0,1)
            tp_bec=capture_bec(tzxy_mf,tp_bec0(kk,:),tr_bec0,0);
            p_bec{jj}(ii,:,kk)=tp_bec;
        end
        tp_halo=mean(p_bec{jj}(ii,:,:),3);
        p_halo{jj}(ii,:)=tp_halo;
        zxy0{ii,jj}=tzxy_mf-tp_halo;
        p_bec0{jj}(ii,:,:)=p_bec{jj}(ii,:,:)-tp_halo;
    end
end

% % DEBUG
% figure(h_zxy_raw);
% hold on;
% plot_zxy(p_halo,[],50,'mk');

h_zxy0=figure();
plot_zxy(zxy0,3e5);
hold on;
for ii=1:2
    for jj=1:2
        plot_zxy({p_bec0{ii}(:,:,jj)},[],100,'k');
    end
end
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
view([0,0]);



%% filter data
zxy0_filt=zxy0;     % initialise filtered data

%%% BEC + thermal
r_thermal=configs.filt.r_ball;
for ii=1:n_shot
    for jj=1:n_mf
        for kk=1:2      % pair of marker BECs
            tp_bec0=p_bec0{jj}(ii,:,kk);
            % get atoms outside ball centered at BEC
            zxy0_filt{ii,jj}=cropBall(zxy0_filt{ii,jj},r_thermal,tp_bec0,false);
        end
    end
end

% DEBUG
h_zxy0_filt_1=figure();
plot_zxy(zxy0_filt,3e5);
hold on;
for ii=1:n_mf
    for jj=1:2
        plot_zxy({p_bec0{ii}(:,:,jj)},[],100,'k');
    end
end
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
view([0,0]);


%%% radial 
% here we rely on "average radius" of halo determined by the mean marker
% BEC locations

% estimate halo radius from marker BECs
r_crop=configs.filt.r_crop;
r_halo=NaN(n_shot,n_mf);
for ii=1:n_mf
    r_halo(:,ii)=0.5*vnorm(p_bec0{ii}(:,:,2)-p_bec0{ii}(:,:,1));
end
r_halo_avg=mean(r_halo,1);

% filter
for ii=1:n_mf
    tr_lim=r_crop*r_halo_avg(ii);       % absolute norm limits
	zxy0_filt(:,ii)=cfilter_norm(zxy0_filt(:,ii),tr_lim(1),tr_lim(2));
end

% DEBUG
h_zxy0_filt_2=figure();
plot_zxy(zxy0_filt,3e5);
hold on;
for ii=1:n_mf
    for jj=1:2
        plot_zxy({p_bec0{ii}(:,:,jj)},[],100,'k');
    end
end
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
view([0,0]);


%%% polar caps
% build box window for halo between caps
z_cap=configs.filt.z_cap;
window_z_filt=cell(1,n_mf);
for ii=1:n_mf
    window_z_filt{ii}={z_cap*r_halo_avg(ii)*[-1,1],[],[]};
end

% filter
for ii=1:n_mf
    tbox_lim=window_z_filt{ii};
    zxy0_filt(:,ii)=cellfun(@(x) boxcull(x,tbox_lim),zxy0_filt(:,ii),'UniformOutput',false);
end

% DEBUG
h_zxy0_filt_3=figure();
plot_zxy(zxy0_filt,3e5);
hold on;
for ii=1:n_mf
    for jj=1:2
        plot_zxy({p_bec0{ii}(:,:,jj)},[],100,'k');
    end
end
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
view([0,0]);



%% k-space and distortion cancellation
k_halo=zxy0_filt;       % initialise atoms in k-space (i.e. atoms lie on unit-sphere)
v_ellip=cell(n_mf,1);

% ellipsoid fit to sphere
for ii=1:n_mf
    [k_halo(:,ii),v_ellip{ii}]=map2usph(k_halo(:,ii));
end
v_ellip=[v_ellip{:}];   % form as struct array
    

% VIS
scatter_halo(k_halo);


%% transform raw spatial clouds
k_all=cell(size(zxy0));
for ii = 1:n_mf
    tV = v_ellip(ii);       
    k_all(:,ii) = cellfun(@(x) ellip2usph(x,tV.cent,tV.rad,tV.vaxis),zxy0(:,ii),...
        'UniformOutput',false);
end

% % VIS
% scatter_halo(k_all);


%% filter
% radial (in/out)
k_all_filt = cfilter_norm(k_all,0.5,1.5);

% hotspot
det_cent=[0,0,0];
hotspot_rad=1.05;
k_all_filt = cellfun(@(k) cylindercull(k,det_cent,[hotspot_rad,Inf],1),k_all_filt,'uni',0);


%%% VIS
scatter_halo(k_all_filt(1:1000,:));



%% categorise data by exp-params
k_all_par=cell(1,nparam);
if nparam>1
    for ii=1:nparam
        k_all_par{ii}=k_all_filt(b_paramset(:,ii),:);      % get all halo data
        %from this param-set and store
    end
else
    % nparam is 1 --> no scan. collate all shots
    k_all_par{1}=k_all_filt;
end

% DEBUG
figure;
for ii=1:nparam
    subplot(1,nparam,ii);
    plot_zxy(k_all_par{ii});
    
    axis equal;
    xlabel('kx');
    ylabel('ky');
    zlabel('kz');
    view([0,0]);
end

% to XYZ coord
k_xyz_all_par = cellfun(@(C) cellfun(@(zxy) zxy2xyz(zxy),C,'uni',0),k_all_par,'uni',0);


%% PUBLICATION: z-measurement (theta=0) result
%%% get data
idx_vis = 1;        % corresponds to theta=0
K = cell_vertcat(k_xyz_all_par{idx_vis});     % REDUCED DATA to xyz


%%% pretty plot
% config
vis_config.FigUnits = 'centimeters';
vis_config.FigPos = [0,0,10,10];
vis_config.FigColor = 'k';

vis_config.Marker = 'o';
vis_config.SizeData = 5;
vis_config.MarkerFaceAlpha = 0.33;
vis_config.MarkerFaceColor = {[0,0,1], [1,0,0]};

% PLOT
H = figure('Name','spinmom_reconstruction');
H.Units = vis_config.FigUnits;
H.Position = vis_config.FigPos;

H.Renderer = 'opengl';  % 'painters'
% H.Renderer = 'painters';
H.InvertHardcopy = 'off';

% set(H,'Color',vis_config.FigColor);       % comment out to make bgd white
% P = scatter3(xyz_filt(:,1),xyz_filt(:,2),xyz_filt(:,3));

idx_vis = 1:10000;

for ii=1:2
    txyz = K{ii};

    P = scatter3(txyz(idx_vis,1),txyz(idx_vis,2),txyz(idx_vis,3));
%     P = scatter3(txyz(:,1),txyz(:,2),txyz(:,3));
    
    set(P,'Marker',vis_config.Marker);
    set(P,'SizeData',vis_config.SizeData);
    set(P,'MarkerEdgeColor','none');
    set(P,'MarkerFaceColor',vis_config.MarkerFaceColor{ii});
    set(P,'MarkerFaceAlpha',vis_config.MarkerFaceAlpha);
    
    hold on;
end

axis equal;
axis off;


% NOTE: make it pretty and save by print_pnghr



% %%% camera
% cpos = [5.1133  -29.5230   15.0022];
% campos(cpos);

