%% Magnetic metrology: two-pulse with delay (Ramsey)
%   approx optimal polarization to eliminate sigma+ component
%   global rotation - beam very large compared to spatial distribution of atoms
%   rotation on (1,1) halo source
%
%   
%
%   NOTE:
%       * significant instability in PUSH/NULLER-SO sequence: mf=1 position variability
%
% 20180323 first
% 20181127 revisited
% DKS
%

%% CONFIGS
% raw data ------------------------------------------------------
run('config_v1');

% Phys-constants ------------------------------------------------
Cphys=physConsts;
C_gymag=Cphys.He_gymag;         % He* gyromagnetic ratio (gamma) [Hz/G]

% vis -----------------------------------------------------------
% publication
config_fig = loadFigureConfig;      % load template

config_fig.pf_alpha=0.33;       % patch face alpha
config_fig.mark_typ={'^','+','o','d'};

config_fig.col_theme=viridis(3);
config_fig.coll_theme=colshades(config_fig.col_theme);

% misc
f_units='normalized';
f_pos=[0.2,0.2,0.2,0.3];
f_pos_wide=[0.2,0.2,0.25,0.3];

fontsize=12;
ax_lwidth=1.2;

%% Get experimental params
% load logfile
param_log=load_logfile(configs.path.paramlog);
param_array = paramlog2array(param_log);
%   format: ID, T_PULSE, PHI_DELAY, T_DELAY

% get unique param-vecs and tag each shot with param-ID
[params,~,Ipar] = unique(flip(param_array(:,2:end),2),'rows');
params=flip(params,2);
% sorted such that the Ith param (vec) can be binned into
%   T_PULSE X PHI_DELAY X T_DELAY (3D) shaped array

param_id=param_array(:,1);
nparam=size(params,1);      % number of unique param-set

par_comp=arrayfun(@(c) unique(params(:,c)), 1:size(params,2),'UniformOutput',false);    % unique param-vec components
ncomp=cellfun(@(x) numel(x), par_comp);  % unique num vals for each comps in the set of param vecs
npar_pulse=ncomp(1);
npar_phidelay=ncomp(2);
npar_tdelay=ncomp(3);

% group shot-ids by exp-param
id_in_param=cell(ncomp);
for ii=1:nparam
    id_in_param{ii}=param_id(Ipar==ii);
end


%% load txy
[txy,fout]=load_txy(configs.load.path,configs.load.id,configs.load.window,configs.load.mincount,configs.load.maxcount,[],1,0,0);

zxy=txy2zxy(txy,vz);
shot_id=fout.id_ok;

% build boolean array selector for scanned param-set
b_paramset=cellfun(@(idx) ismember(shot_id,idx),...
    id_in_param,'UniformOutput',false);
% b_paramset=horzcat(b_paramset{:});


h_zxy_raw=figure;
plot_zxy(zxy,1e5);
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
view([0,0]);

%% VIS: 3D scatter raw data -----------------------------------
idx_Ramsey=1;       % pi/2 pulse
idx_phi0=1;         % phase delay = 0
idx_tau_vis = 1:4:16;   % indices of taus to visualise
n_tau_vis=length(idx_tau_vis);
b_sel=b_paramset(idx_Ramsey,idx_phi0,:);


% ZXY with Z-values relative to initial v_z=0 free-fall
t_trapSO=0;         % time at trap switchoff
t_tof=0.416;        % free-fall TOF till detection
zxy_zrel=cellfun(@(r) r - [vz*(t_trapSO+t_tof),0,0],zxy,'uni',0);
% TODO code TXY --> ZXY where rel-Z is already calculated
%   relZ will be proportional to v_z at trapSO

for ii = 1:n_tau_vis
    H_zxy_raw_3d = figure('Name','zxy_raw_3d','Units',config_fig.units,'Position',[0 0 8.6 12],'Renderer',config_fig.rend);
    hold on;
    
    % select Ramsey
    idx_tau=idx_tau_vis(ii);
%     b_sel=b_paramset{idx_Ramsey,idx_phi0,idx_tau};
    zxy_zrel_sel=zxy_zrel(b_sel{idx_tau});
    
    s = vis_zxy_3d(1e3*zxy_zrel_sel{1});
    
    ax=gca;
    set(ax,'XLim',1e3*configs.load.window{2});
    set(ax,'YLim',1e3*configs.load.window{3});
    xlabel('$x$ (mm)');
    ylabel('$y$ (mm)');
    zlabel('$z^*$ (mm)');
    
    titlestr = sprintf('%s = %0.3g %s','$\tau$',1e6*par_comp{3}(idx_tau),'$\mu$s');
    title(titlestr);
    
    H_zxy_raw_3d.Name=strcat(H_zxy_raw_3d.Name,'_',num2str(1e6*par_comp{3}(idx_tau)));
end

%% VIS: projection raw data -----------------------------------
img_ax = 'xz';
img_lim = {configs.load.window{2}, 1e-3*[-200,100]};
img_pitch = {1e-3, 1e-3};

I=cell(n_tau_vis,1);

% side-by-side ------------------------------------------
H = figure;
for ii = 1:n_tau_vis
    ax=subplot(1,n_tau_vis,ii);
    
    % select data subset (Ramsey) @ particular tau
    idx_tau=idx_tau_vis(ii);
%     b_sel=b_paramset{idx_Ramsey,idx_phi0,idx_tau};
    zxy_zrel_sel=zxy_zrel(b_sel{idx_tau});
    
    ZXY_this=cat(1,zxy_zrel_sel{:});
    [I{ii},xx,yy] = zxy2img(ZXY_this,img_ax,img_lim,img_pitch);
    imagesc(ax,'XData',xx,'YData',yy,'CData',log10(I{ii}));
    
    axis tight;
    axis equal;
    axis off;
    
    if ii==length(idx_tau_vis)
        pos0=ax.Position;   %original axes position
        % colorbar
        cbar=colorbar('eastoutside');
        cbar.Label.String='log_{10}(atom number)';
        %return axes pos
        ax.Position=pos0;
    end
end

% juxtaposed images ------------------------------------------
Icat = cat(2,I{:});     % along horizontal

H_xz_raw_proj = figure('Name','xz_raw_proj','Units',config_fig.units,'Position',[0 0 13 12],'Renderer',config_fig.rend);
imagesc('XData',xx*n_tau_vis,'YData',yy,'CData',log10(Icat));

ax=gca;
axis tight;
axis equal;
axis off;
pos0=ax.Position;   %original axes position

cbar=colorbar('eastoutside');       % colorbar
cbar.Label.String='log_{10}(atom number)';
ax.Position=pos0;   %return axes pos
boxpos0=plotboxpos(ax);
cbar.Position(3)=cbar.Position(3)/3;
cbar.Position(1)=boxpos0(1)+boxpos0(3)+cbar.Position(3);
colormap('inferno');

% annotation
for ii=1:n_tau_vis
    idx_tau=idx_tau_vis(ii);
    tt = 1e6*par_comp{3}(idx_tau);
    tx = xx(1)*n_tau_vis + (ii-1+0.5)*diff(headtail(xx));
    ty = yy(end) + 0.02*diff(headtail(yy));
    
    textstr = sprintf('%s = %0.3g %s','$\tau$',tt,'$\mu$s');
    text(tx,ty,textstr,'HorizontalAlignment','center','VerticalAlignment','bottom');
end


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

% DEBUG
% TODO - pretty this figure for publication
figure(h_zxy_raw);
hold on;
plot_zxy(p_halo,[],50,'mkg');

h_zxy0=figure();
plot_zxy(zxy0,3e5);
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


%% VIS
% plot spin (mJ) distinguished result + halo centered (oscillations
% cancelled)

% 3D scatter -----------------------------------
for ii=1:2  %mJ    
    H=figure;
    
    for jj=1:n_tau_vis
        subplot(1,n_tau_vis,jj);
        
        idx_tau=idx_tau_vis(jj);
        
        Z = zxy0{find(b_sel{idx_tau},1),ii};
        vis_zxy_3d(1e3*Z);
        
        ax=gca;
        
        ax_lim=[-35,35];
        set(ax,'XLim',ax_lim);
        set(ax,'YLim',ax_lim);
        set(ax,'ZLim',ax_lim);
    end
end

%% VIS: Images
% 2D projection  ---------------------------------------------------
b_cat = any(cat(2,b_sel{:}),2);     % any shot belonging to this experiment
img_lim = {35e-3*[-1,1],35e-3*[-1,1]};
img_pitch = {0.2e-3,0.2e-3};

% XY
H=figure('Name','mJ_xy_proj','Units',config_fig.units,'Position',[0 0 6 12],'Renderer',config_fig.rend);
for ii=1:2
    subplot(2,1,ii);
    
    Z=cat(1,zxy0{b_cat,ii});        % concatenate all mJ=ii atoms
    
    [I,xx,yy]=zxy2img(Z,'xy',img_lim,img_pitch);
    
    imagesc('XData',1e3*xx,'YData',1e3*yy,'CData',log10(I));
    
    ax=gca;
    axis tight;
    axis equal;
    colormap('inferno');
    xlabel('$x$ (mm)');
    ylabel('$y$ (mm)');
    
%     % colorbar
%     pos0=ax.Position;   %original axes position
%     cbar=colorbar('eastoutside');
%     cbar.Label.String='log_{10}(atom number)';
%     ax.Position=pos0;   %return axes pos
%     boxpos0=plotboxpos(ax);
%     cbar.Position(3)=cbar.Position(3)/3;
%     cbar.Position(1)=boxpos0(1)+boxpos0(3)+cbar.Position(3);
%     colormap('inferno');
end

% XZ
H=figure('Name','mJ_xz_proj','Units',config_fig.units,'Position',[0 0 6 12],'Renderer',config_fig.rend);
for ii=1:2
    subplot(2,1,ii);
    
    Z=cat(1,zxy0{b_cat,ii});        % concatenate all mJ=ii atoms
    
    [I,xx,yy]=zxy2img(Z,'xz',img_lim,img_pitch);
    
    imagesc('XData',1e3*xx,'YData',1e3*yy,'CData',log10(I));
    
    ax=gca;
    axis tight;
    axis equal;
    colormap('inferno');
    xlabel('$x$ (mm)');
    ylabel('$z$ (mm)');
    
%     % colorbar
%     pos0=ax.Position;   %original axes position
%     cbar=colorbar('eastoutside');
%     cbar.Label.String='log_{10}(atom number)';
%     ax.Position=pos0;   %return axes pos
%     boxpos0=plotboxpos(ax);
%     cbar.Position(3)=cbar.Position(3)/3;
%     cbar.Position(1)=boxpos0(1)+boxpos0(3)+cbar.Position(3);
%     colormap('inferno');
end

%% XY-slices -----------------------------------------
img_lim = {35e-3*[-1,1],35e-3*[-1,1]};
img_pitch = {0.4e-3,0.4e-3};

z0 = linspace(-26e-3,26e-3,5);
dz = 52e-3/20;

H=figure('Name','mJ_xy_slice','Units',config_fig.units,'Position',[0 0 26 12],'Renderer',config_fig.rend);

lim_zslice = arrayfun(@(z) z + dz*[-1,1],z0,'uni',0);

for ii=1:2
    Z=cat(1,zxy0{b_cat,ii});        % concatenate all mJ=ii atoms
    
    Zsliced = cellfun(@(z) boxcull(Z,{z,[],[]}),lim_zslice,'uni',0);
    
    I=cellfun(@(z) zxy2img(z,'xy',img_lim,img_pitch),Zsliced,'uni',0);
    
    Icat=cat(2,I{:});      % stitch images horizontally
    
    ax=subplot(2,1,ii);
    imagesc('XData',length(z0)*img_lim{1},...
        'YData',img_lim{2},...
        'CData',log10(Icat));
    
    
    axis equal;
    axis tight;
    axis off;
%     ax.XTickLabel=[];
%     ax.YTickLabel=[];
    colormap('inferno');
    
    % colorbar
    pos0=ax.Position;   %original axes position
    cbar=colorbar('eastoutside');       
    cbar.Label.String='log_{10}(atom number)';
    ax.Position=pos0;   %return axes pos
    boxpos0=plotboxpos(ax);
    cbar.Position(3)=cbar.Position(3)/3;
    cbar.Position(1)=boxpos0(1)+boxpos0(3)+cbar.Position(3);
    colormap('inferno');
    
    % annotation
    if ii==1
        for jj=1:length(z0)
            xx=img_lim{1};
            yy=img_lim{2};
            
            tz = 1e3*z0(jj);
            tx = xx(1)*length(z0) + (jj-1+0.5)*diff(headtail(xx));
            ty = yy(end) + 0.02*diff(headtail(yy));
            
            textstr = sprintf('$z$ = %0.3g mm',tz);
            text(tx,ty,textstr,'HorizontalAlignment','center','VerticalAlignment','bottom');
        end
    end
end


%% filter data
zxy0_filt=zxy0;     % initialise filtered data

% non-halo stuff
r_bec_capt=configs.crop.bec;
zxy0_bec=cell(n_shot,n_mf,2);


%%% BEC + thermal
r_thermal=configs.filt.r_ball;
for ii=1:n_shot
    for jj=1:n_mf
        for kk=1:2      % pair of marker BECs
            tp_bec0=p_bec0{jj}(ii,:,kk);
            % get atoms outside ball centered at BEC
            zxy0_filt{ii,jj}=cropBall(zxy0_filt{ii,jj},r_thermal,tp_bec0,false);
            
            % get BEC (relatively tight)
            zxy0_bec{ii,jj,kk}=cropBall(zxy0{ii,jj},r_bec_capt,tp_bec0,true);     % BEC
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


%% VIS - ellipsoid fit (3D)
for ii=1:2
    Z = zxy0_filt(b_cat,ii);        % truncated halo raw data points
    tV=v_ellip(ii);      % ellipsoid data
    
    H = vis_ellipsoid_fit(zxy2xyz(cat(1,Z{:})),tV.cent,tV.rad,tV.vaxis,tV.v);
    
    figname=strcat('ellipsoid_fit_',num2str(ii));
    H.Name=figname;
    
    % print this via HR PNG printer
    %   print(gcf,'foo.png','-dpng','-r600');
    % or
    %   print_pnghr(gcf);
end


%% VIS - ellipsoid fit (2D slice)
% configs ------------------------------
img_lim = {35e-3*[-1,1],35e-3*[-1,1]};
img_pitch = {0.4e-3,0.4e-3};

zf = 0.75;       % factor from max radii
z0 = linspace(-26e-3*zf,26e-3*zf,5);
dz = 52e-3/20;

lim_zslice = arrayfun(@(z) z + dz*[-1,1],z0,'uni',0);

ax_ngrid=100;
[xx,yy]=ndgrid(linspace(img_lim{1}(1),img_lim{1}(2),ax_ngrid),...
    linspace(img_lim{2}(1),img_lim{2}(2),ax_ngrid));

H=figure('Name','ellipsoid_fit_contour','Units',config_fig.units,'Position',[0 0 26 9.5],'Renderer',config_fig.rend);
hax=tight_subplot(2,length(z0),0,0.05,0.1);
counter=1;
for ii=1:2
    Z = cat(1,zxy0_filt{b_cat,ii});         % concatenate all mJ=ii atoms (truncated)
    % loop through Z-slices
    Zsliced = cellfun(@(z) boxcull(Z,{z,[],[]}),lim_zslice,'uni',0);
    % slice density image (XY)
    I=cellfun(@(z) zxy2img(z,'xy',img_lim,img_pitch),Zsliced,'uni',0);
    
    for jj=1:length(z0)
        ax=hax(counter);
        axes(ax);
        
        hold on;
        
        imagesc('XData',img_lim{1},'YData',img_lim{2},'CData',log10(I{jj}));
        colormap('inferno')
        
        % fitted ellipsoid contour
        f_ellip=ellipsoidEqn(v_ellip(ii).v,xx,yy,z0(jj)*ones(size(xx)));
        [c,h]=contour(xx,yy,f_ellip,'LevelList',[0]);   % contour at this slice
        h.LineColor='g';
        h.LineStyle='--';
        
        axis equal;
        axis tight;
        axis off;
        
        %title
        if ii==1
            textstr = sprintf('$z$ = %0.3g mm',1e3*z0(jj));
            title(textstr);
        end
        
        % colorbar
        if jj==length(z0)
            pos0=ax.Position;   %original axes position
            cbar=colorbar('eastoutside');
            cbar.Label.String='log_{10}(atom number)';
            ax.Position=pos0;   %return axes pos
            boxpos0=plotboxpos(ax);
            cbar.Position(3)=cbar.Position(3)/3;
            cbar.Position(1)=boxpos0(1)+boxpos0(3)+cbar.Position(3);
            colormap('inferno');
        end
        
        counter=counter+1;
    end
end


%% filter post-processed data
%%% radial
r_crop2=configs.filt2.r_crop;
k_halo_filt=cfilter_norm(k_halo,r_crop2(1),r_crop2(2));

%%% z-cap
% build box window for halo between caps
z_cap2=configs.filt2.z_cap;
window_z_filt2=cell(1,n_mf);
for ii=1:n_mf
    window_z_filt2{ii}={z_cap2*[-1,1],[],[]};
end

% filter
for ii=1:n_mf
    tbox_lim=window_z_filt2{ii};
    k_halo_filt(:,ii)=cellfun(@(x) boxcull(x,tbox_lim),k_halo_filt(:,ii),'UniformOutput',false);
end


% DEBUG
scatter_halo(k_halo_filt);


%% transform raw spatial clouds
k_all=cell(size(zxy0));
for ii = 1:n_mf
    tV = v_ellip(ii);       
    k_all(:,ii) = cellfun(@(x) ellip2usph(x,tV.cent,tV.rad,tV.vaxis),zxy0(:,ii),...
        'UniformOutput',false);
end

% VIS
scatter_halo(k_all);

%%% filter radially
k_all_filt = cfilter_norm(k_all,configs.filt2.r_crop(1),configs.filt2.r_crop(2));

% VIS
scatter_halo(k_all_filt);


%% select data to analyse
% if ~isfield(configs,'roi')
%     warning('Region to analyse has not been specified. Default: scattering halo "k_halo_filt".');
%     k_roi=k_halo_filt;
% elseif strcmp(configs.roi,'halo')
%     k_roi=k_halo_filt;
% elseif strcmp(configs.roi,'all')
%     k_roi=k_all_filt;
% else
%     error('Unrecognised region of interest to analyse.');
% end


%% categorise data by exp-params
k_halo_par=cell(ncomp);     % halo only
k_all_par=cell(ncomp);      % spherical shell (w/ BEC)
x_bec_par=cell(ncomp);      % BEC counts. spatial coord [m]

for ii=1:nparam
    k_all_par{ii}=k_all_filt(b_paramset{ii},:);   
    k_halo_par{ii}=k_halo_filt(b_paramset{ii},:); 
    
    x_bec_par{ii}=zxy0_bec(b_paramset{ii},:,:);
end


% NOTE: generates a lot of subplots!
% % DEBUG
% figure;
% for ii=1:nparam
%     subplot(1,nparam,ii);
%     plot_zxy(k_par{ii});
%     
%     axis equal;
%     xlabel('kx');
%     ylabel('ky');
%     zlabel('kz');
%     view([0,0]);
% end


%%% END OF PREPROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% clean workspace
clear txy zxy zxy0 zxy0_filt zxy0_bec tzxy tzxy_mf tp_bec tp_bec0 tp_halo tr_bec0 tr_lim;
clear p_bec p_bec0 p_halo r_bec_capt r_crop* r_halo* r_thermal tbox_lim window_z_filt*;
clear h_zxy*;       % clear figs


%% transform to RH coords
% SAVE ORIGINAL 
if ~exist('k_all_orig','var')
    k_all_orig=k_all_par;
    k_halo_orig=k_halo_par;
end

% EXP-coord sys (z against g)
k_all_par=cellfun(@(C) cellfun(@(x) tzxy2RHtzxy2(x),C,'UniformOutput',false),...
    k_all_orig,'UniformOutput',false);      
k_halo_par=cellfun(@(C) cellfun(@(x) tzxy2RHtzxy2(x),C,'UniformOutput',false),...
    k_halo_orig,'UniformOutput',false);      


%% ANALYSIS
tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRE-ANALYSIS: select subset of DATA relevant

tau=par_comp{3};              % T_delay between two pulses [s]
phi=par_comp{2};            % phase delay of 2nd pi/2 pulse (rad)

% Pre-processing
idx_ramsey_exp=1;           % 1: Ramsey (pi/2), 2: pi-pulses
idx_phi0=find(phi==0);      % two-pi/2-pulse (Ramsey) dataset with zero phase delay phi=0

% Ramsey sub-dataset (phi=0)
k_all_ramsey_phi0=squeeze(k_all_par(idx_ramsey_exp,idx_phi0,:));     
k_halo_ramsey_phi0=squeeze(k_halo_par(idx_ramsey_exp,idx_phi0,:));

x_ramsey_phi0=squeeze(x_bec_par(idx_ramsey_exp,idx_phi0,:));    % BEC sub-dataset


%% General halo/collision stuff
% Ramsey (3ms)
d_sep = 3e-3 * 120e-3;      % 3ms expansion at 120 mm/s
sig_psf_r = 35e-6;      % rms width of point-spread function ~ 35um (BEC) (see supplementary in paper)

sig_psf_beta=sig_psf_r/(d_sep/2);   % PSF rms width in angle (rad)
% sig_filt_beta=sig_psf_beta/2;       % filter rms width (rad)

% gradiometry (~1.5 ms)
% sig_psf_gradiometry=sig_psf_beta*2;     % ~twice uncertainty
% 1.7 ms
sig_psf_gradiometry=0.156*pi;     % ~twice uncertainty



%% atom number and population
n_shot=cellfun(@(x) size(x,1),k_all_ramsey_phi0);       % number of shots per tau

% raw #spins: #atom in mJ halos
Nm_all=cellfun(@(x) shotSize(x),k_all_ramsey_phi0,'UniformOutput',false);      
Nm_halo=cellfun(@(x) shotSize(x),k_halo_ramsey_phi0,'UniformOutput',false);      

% #spin in halo
Nm_halo_avg=cellfun(@(n) mean(n,1),Nm_halo,'UniformOutput',false);      % avg num in mJ halo
Nm_halo_avg=vertcat(Nm_halo_avg{:});        % form as array
Nm_halo_std=cellfun(@(n) std(n,0,1),Nm_halo,'UniformOutput',false);     % std 
Nm_halo_std=vertcat(Nm_halo_std{:});
Nm_halo_se=Nm_halo_std./sqrt(n_shot);       % std-err

Ninv_halo=cellfun(@(x) -diff(x,1,2),Nm_halo,'UniformOutput',false);    % "collective spin"/inversion

% tot #atom in halo per tau
N_halo=cellfun(@(n) sum(n,2),Nm_halo,'UniformOutput',false);    % total num det in halo
N_halo_avg=cellfun(@(x) mean(x),N_halo);    % avg total atom num det'd (T avg'd)
N_halo_std=cellfun(@(x) std(x),N_halo);
N_halo_se=N_halo_std./sqrt(n_shot);

% tot #atom in halo averaged over tau
N_halo_avg_exp=mean(vertcat(N_halo{:}));    
N_halo_std_exp=std(vertcat(N_halo{:}));
N_halo_se_exp=N_halo_std_exp/sqrt(sum(n_shot));

%%% population
% spin-components
pm_halo=cellfun(@(n,N) n./N,Nm_halo,N_halo,'UniformOutput',false);

pm_halo_avg=cellfun(@(x) mean(x,1),pm_halo,'UniformOutput',false);
pm_halo_avg=vertcat(pm_halo_avg{:});
pm_halo_std=cellfun(@(x) std(x,0,1),pm_halo,'UniformOutput',false);
pm_halo_std=vertcat(pm_halo_std{:});
pm_halo_se=pm_halo_std./sqrt(n_shot);

% P: population inversion: P = p_up - p_down
P_halo=cellfun(@(n,N) n./N,Ninv_halo,N_halo,'UniformOutput',false);   

P_halo_avg=cellfun(@(x) mean(x),P_halo);
P_halo_std=cellfun(@(x) std(x),P_halo);
P_halo_se=P_halo_std./sqrt(n_shot);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BEC
Nm_bec=cellfun(@(x) shotSize(x),x_ramsey_phi0,'UniformOutput',false);

Nm_bec_avg=cellfun(@(n) mean(n,1),Nm_bec,'UniformOutput',false);
Nm_bec_avg=vertcat(Nm_bec_avg{:});
Nm_bec_std=cellfun(@(n) std(n,0,1),Nm_bec,'UniformOutput',false);     % std 
Nm_bec_std=vertcat(Nm_bec_std{:});
Nm_bec_se=Nm_bec_std./sqrt(n_shot);       % std-err
Ninv_bec=cellfun(@(x) -squeeze(diff(x,1,2)),Nm_bec,'UniformOutput',false);    % "collective spin"/inversion

N_bec=cellfun(@(n) squeeze(sum(n,2)),Nm_bec,'UniformOutput',false);
N_bec_avg=cellfun(@(x) mean(x,1),N_bec,'UniformOutput',false);    % avg total atom num det'd (T avg'd)
N_bec_avg=vertcat(N_bec_avg{:});
N_bec_std=cellfun(@(x) std(x,0,1),N_bec,'UniformOutput',false);
N_bec_std=vertcat(N_bec_std{:});
N_bec_se=N_bec_std./sqrt(n_shot);

P_bec=cellfun(@(n,N) n./N,Ninv_bec,N_bec,'UniformOutput',false);   
P_bec_avg=cellfun(@(x) mean(x),P_bec,'UniformOutput',false);
P_bec_avg=vertcat(P_bec_avg{:});
P_bec_std=cellfun(@(x) std(x),P_bec,'UniformOutput',false);
P_bec_std=vertcat(P_bec_std{:});
P_bec_se=P_bec_std./sqrt(n_shot);


%% VIS: HALO COLLECTIVE RAMSEY SIGNAL
% configs
col_spin={'b','r'};     % mJ = 1, 0

% PLOT =====================================================
figname='halo_ramsey_signal';
h=figure('Name',figname,'Units',config_fig.units,'Position',[0,0,config_fig.pagewidth,6],...
    'Renderer',config_fig.rend);
hold on;

% detected atom number ------------------------------------------------
subplot(1,2,1);
hold on;
for ii=1:2    
    tp_halo=ploterr(1e6*tau,Nm_halo_avg(:,ii),[],Nm_halo_se(:,ii),'o-');
    set(tp_halo(1),'Color',col_spin{ii},'MarkerFaceColor',col_spin{ii},...
        'MarkerSize',config_fig.mark_siz,'LineWidth',config_fig.line_wid);
    set(tp_halo(2),'Color',col_spin{ii},'LineWidth',config_fig.line_wid);
end

ax=gca;
box on;
set(ax,'Layer','Top');
set(ax,'FontSize',config_fig.ax_fontsize);
set(ax,'LineWidth',config_fig.ax_lwid);

axis_snug(ax,[0.05,0.1]);

xlabel('pulse delay $\tau~(\mu s)$');
ylabel('\# detected $N_i$');

% Polarisation ------------------------------------------------
subplot(1,2,2);
hold on;

tp_halo=ploterr(1e6*tau,P_halo_avg,[],P_halo_se,'o-');
set(tp_halo(1),'Color','k','MarkerFaceColor','k',...
    'MarkerSize',config_fig.mark_siz,'LineWidth',config_fig.line_wid);
set(tp_halo(2),'Color','k',...
    'MarkerSize',config_fig.mark_siz,'LineWidth',config_fig.line_wid);

ax=gca;
box on;
set(ax,'Layer','Top');
set(ax,'FontSize',config_fig.ax_fontsize);
set(ax,'LineWidth',config_fig.ax_lwid);

xlabel('pulse delay $\tau~(\mu s)$');
ylabel('polarisation');
axis_snug(ax,[0.05,0.1]);

%% VIS: BEC COLLECTIVE RAMSEY SIGNAL
% configs
i_bec=1:2;
mrk_bec={'v','^'};      % marker types for BEC: triangle DOWN/UP
col_spin={'b','r'};     % mJ = 1, 0

% PLOT =====================================================
figname='bec_ramsey_signal';
h=figure('Name',figname,'Units',config_fig.units,'Position',[0,0,config_fig.pagewidth,6],...
    'Renderer',config_fig.rend);
hold on;

% detected atom number ----------------------------------------
subplot(1,2,1);
hold on;

for ii=1:2      % SPIN component
    for jj=i_bec    % BEC North/South
        tp_bec=ploterr(1e6*tau,Nm_bec_avg(:,ii,jj),[],Nm_bec_se(:,ii,jj),'--');
        set(tp_bec(1),'Color',col_spin{ii},'MarkerFaceColor','w','Marker',mrk_bec{jj},...
            'MarkerSize',config_fig.mark_siz,'LineWidth',config_fig.line_wid);
        set(tp_bec(2),'Color',col_spin{ii},'LineWidth',config_fig.line_wid);
        
        if jj==1
            tp_bec(1).MarkerFaceColor=tp_bec(1).Color;
        end
    end
end

ax=gca;
box on;
set(ax,'Layer','Top');
set(ax,'FontSize',config_fig.ax_fontsize);
set(ax,'LineWidth',config_fig.ax_lwid);

axis_snug(ax,[0.05,0.1]);

xlabel('pulse delay $\tau~(\mu s)$');
ylabel('\# detected $N_i$');


% Polarisation -----------------------------------------
subplot(1,2,2);
hold on;
% bec
for jj=i_bec
    tp_bec=ploterr(1e6*tau,P_bec_avg(:,jj),[],P_bec_se(:,jj),'--');
    set(tp_bec(1),'Color','k','MarkerFaceColor','w','Marker',mrk_bec{jj},...
        'MarkerSize',config_fig.mark_siz,'LineWidth',config_fig.line_wid);
    set(tp_bec(2),'Color','k',...
        'MarkerSize',config_fig.mark_siz,'LineWidth',config_fig.line_wid);
    if jj==1
        tp_bec(1).MarkerFaceColor=tp_bec(1).Color;
    end
end

ax=gca;
box on;
set(ax,'Layer','Top');
set(ax,'FontSize',config_fig.ax_fontsize);
set(ax,'LineWidth',config_fig.ax_lwid);

xlabel('pulse delay $\tau~(\mu s)$');
ylabel('polarisation');
axis_snug(ax,[0.05,0.1]);


%% Ramsey model
% polarisation/population inversion ---------------------------------------
Pramsey_mdl='y~amp*cos(om*x+phi)';
Pramsey_cname={'amp','om','phi'};
Pramsey_par0=[1,2*pi*1.5,pi];

idx_om_Pramsey=2;     %find(cellfun(@(s) strcmp(s,'om'),Pramsey_cname)==1);

% domain x to evaluate fit
tt_fit=1e6*linspace(min(tau),max(tau),1e3);     % x-axis range for fitted curve

%% collate all data x-values
tt=arrayfun(@(x,n) x*ones(n,1),tau,n_shot,'UniformOutput',false);
tt=vertcat(tt{:});

%% fit to halo-avgd data
%%% Polarisation --------------------------------------------------------
% collate data to fit
yy=vertcat(P_halo{:});

% fit
Pramsey_fit_halo=fitnlm(1e6*tt,yy,Pramsey_mdl,Pramsey_par0,'CoefficientNames',Pramsey_cname);

% get fit params
Pramsey_fpar_halo=Pramsey_fit_halo.Coefficients.Estimate;
Pramsey_fparerr_halo=Pramsey_fit_halo.Coefficients.SE;

% get freq
om_Pramsey_halo=Pramsey_fpar_halo(idx_om_Pramsey);
omerr_Pramsey_halo=Pramsey_fparerr_halo(idx_om_Pramsey);

% magnetic field
B_Pramsey_halo=1e6*om_Pramsey_halo/(2*pi*C_gymag);
Berr_Pramsey_halo=1e6*omerr_Pramsey_halo/(2*pi*C_gymag);


%% VIS: Ramsey fringe (pop-inv)
figname='Pramsey_fringe_halo';
h=figure('Name',figname,'Units',f_units,'Position',[0.2,0.2,0.5,0.2],'Renderer',config_fig.rend);
hold on;

pexp=ploterr(1e6*tau,P_halo_avg,...
    [],P_halo_se,'o','hhxy',0);
set(pexp(1),'MarkerFaceColor',config_fig.coll_theme(1,:),'MarkerEdgeColor',config_fig.col_theme(1,:));
set(pexp(2),'Color',config_fig.col_theme(1,:));

yy_fit=feval(Pramsey_fit_halo,tt_fit);
pfit=plot(tt_fit,yy_fit,'LineStyle',config_fig.line_sty{1},'Color',config_fig.coll_theme(1,:));
uistack(pfit,'bottom');

titlestr=sprintf('%s %0.3g(%0.1g) MHz','$f_L=$',om_Pramsey_halo/(2*pi),omerr_Pramsey_halo(idx_phi0)/(2*pi));
title(titlestr);

% annotate subplot
ax=gca;
box on;
ax_ylim=ax.YLim;
ax_xlim=ax.XLim;

set(ax,'Layer','Top');
ax.FontSize=fontsize;
ax.LineWidth=ax_lwidth;
xlabel('pulse delay $\tau~(\mu s)$');
ylabel('$P$');

% % save fig
% if do_save_figs
%     savefigname=sprintf('fig_%s_%s',figname,getdatetimestr);
%     fpath=fullfile(dir_save,savefigname);
%     
%     saveas(h,strcat(fpath,'.fig'),'fig');
%     print(h,strcat(fpath,'.svg'),'-dsvg');
% end


%% atom number distribution: nk_m
%%% Spatial zones
% construct spatial zones at latlon grid + solid angle
% half-cone angle of integration bin (rad)
% alpha=sig_psf_beta;          % limiting spatial resolution 
alpha=sig_psf_gradiometry;      % bin-size like gradiometry (rad)
lim_az=[-pi,pi];        % no inversion symmetry
phi_max=pi/2;           
lim_el=[-phi_max,phi_max];

n_az=200;                	% equispaced bins
n_el=100;

% % QUICK DEBUG
% n_az=100;
% n_el=50;


az_disp=deg2rad(-180:90:90);     % azim sections (great circles) to display
% el_disp=deg2rad(-30:30:30);    % elev/lat zones to display

az=linspace(lim_az(1),lim_az(2),n_az+1);
az=az(1:end-1);
el=linspace(lim_el(1),lim_el(2),n_el);

[~,iaz_disp]=arrayfun(@(x) min(abs(az-x)),az_disp);     % idx to displayable azim
naz_disp=length(iaz_disp);

[gaz,gel]=ndgrid(az,el);    % az-el grid
n_zone=numel(gaz);

% idx to ~0 angle (equator for elev)
[~,iaz_0]=min(abs(az));
[~,iel_0]=min(abs(el));


%% n,P distribution
%%% atom number
Nm_k=cell(numel(tau),1);     % #atoms in k,mj-zone categorise by exp param

% evaluate
for ii=1:numel(k_all_ramsey_phi0)
    tk=k_all_ramsey_phi0{ii};        % exp data for this param
    
    % get num in zone/shot
    tN=arrayfun(@(th,ph) cellfun(@(x) size(inCone(x,th,ph,alpha),1),tk),...
        gaz,gel,'UniformOutput',false);     
    
    % format into multi-dim array: AZ X EL X SHOT X MJ
    ttN=NaN(cat(2,size(gaz),size(tN{1})));
    for jj=1:n_zone
        [iaz,iel]=ind2sub([n_az,n_el],jj);
        
        ttN(iaz,iel,:,:)=tN{iaz,iel};
    end
    Nm_k{ii}=ttN;
end

Ninv_k=cellfun(@(x) -diff(x,1,4),Nm_k,'UniformOutput',false);
N_k=cellfun(@(x) sum(x,4),Nm_k,'UniformOutput',false);

%%% population 
% polarisation ------------------------------------------------------------
P_k=cellfun(@(n,N) n./N,Ninv_k,N_k,'UniformOutput',false);
P_k_avg=cellfun(@(x) squeeze(mean(x,3,'omitnan')),P_k,'UniformOutput',false);
P_k_std=cellfun(@(x) squeeze(std(x,0,3,'omitnan')),P_k,'UniformOutput',false);

% multidim array: PHI_DELAY X T_DELAY X AZ X EL
P_k_avg=cell2mat(cellfun(@(a) reshape(a,[1,size(a)]),P_k_avg,'UniformOutput',false));
P_k_std=cell2mat(cellfun(@(a) reshape(a,[1,size(a)]),P_k_std,'UniformOutput',false));
P_k_se=P_k_std./sqrt(n_shot);        


%% equatorial integrated number
Nm_eq=cellfun(@(x) squeeze(x(:,iel_0,:,:)),Nm_k,'UniformOutput',false);     % k in equator
Nm_eq=cellfun(@(x) squeeze(sum(x,1)),Nm_eq,'UniformOutput',false);          % equator-integrated

Ninv_eq=cellfun(@(x) -diff(x,1,2),Nm_eq,'UniformOutput',false);
N_eq=cellfun(@(x) sum(x,2),Nm_eq,'UniformOutput',false);

P_eq=cellfun(@(n,N) n./N,Ninv_eq,N_eq,'UniformOutput',false);
P_eq_avg=cellfun(@(x) mean(x),P_eq);
P_eq_std=cellfun(@(x) std(x),P_eq);
P_eq_se=P_eq_std./sqrt(n_shot);


%% k-resolved fit
%%% pop-inversion --------------------------------------------------------
Pramseyk_fit=cell(n_az,n_el);   % dim: AZ X EL
for ii=1:n_zone
    [iaz,iel]=ind2sub([n_az,n_el],ii);
    
    % collate data to fit
    yy=cellfun(@(x) squeeze(x(iaz,iel,:)),P_k,'UniformOutput',false);
    yy=vertcat(yy{:});
    
    % fit
    Pramseyk_fit{iaz,iel}=fitnlm(1e6*tt,yy,Pramsey_mdl,...
        Pramsey_par0,'CoefficientNames',Pramsey_cname);
end

% get fit params
Pramseyk_fpar=arrayfun(@(I) cellfun(@(f) f.Coefficients.Estimate(I),Pramseyk_fit),...
    1:numel(Pramsey_cname),'UniformOutput',false);
Pramseyk_fpar=cat(ndims(Pramseyk_fpar{1})+1,Pramseyk_fpar{:});

Pramseyk_fparerr=arrayfun(@(I) cellfun(@(f) f.Coefficients.SE(I),Pramseyk_fit),...
    1:numel(Pramsey_cname),'UniformOutput',false);
Pramseyk_fparerr=cat(ndims(Pramseyk_fparerr{1})+1,Pramseyk_fparerr{:});
% dim: AZ X EL X PAR#

% get freq
omk_Pramsey=Pramseyk_fpar(:,:,idx_om_Pramsey);
omerrk_Pramsey=Pramseyk_fparerr(:,:,idx_om_Pramsey);
% dim: AZ X EL

% % outlier to NaN
% b_outlier=isoutlier(omk_Pramsey);    % find outlier by Med Abs Dev
% omk_Pramsey_raw=omk_Pramsey;      % store original raw data set
% omerrk_Pramsey_raw=omerrk_Pramsey;
% omk_Pramsey(b_outlier)=NaN;          % outlier --> NaN
% omerrk_Pramsey(b_outlier)=NaN;

% omk_Pramsey_0=squeeze(mean(omk_Pramsey,1,'omitnan'));   % mode-resolved L-freq avgd thru PHI_DELAY
% n_samp=squeeze(sum(~b_outlier,1));       % num exp samples to average over
% omerrk_Pramsey_0=sqrt(squeeze(sum(omerrk_Pramsey.^2,1,'omitnan')))./n_samp;

%%% magnetic field 
%%%% far-field
Bk_Pramsey=1e6*omk_Pramsey/(2*pi*C_gymag);
Berrk_Pramsey=1e6*omerrk_Pramsey/(2*pi*C_gymag);

% equatorial slice
Bk_eq=Bk_Pramsey(:,iel_0);
Bkerr_eq=Berrk_Pramsey(:,iel_0);      


% %%%% interrogation
% % METHOD: Evaluate magnetic field estimated at interrogation time
% %   Gaussian convolution: see SMs for determination of length scale
% %       ERROR by conv of variance
% % Br=gaussfilt_sph(Bk_Pramsey,gaz,gel,gaz,gel,sig_psf_beta);      
% % Brerr=sqrt(gaussfilt_sph(Berrk_Pramsey.^2,gaz,gel,gaz,gel,sig_psf_beta));
% % 
% Br=gaussfilt_sph(Bk_Pramsey,gaz,gel,gaz,gel,sig_filt_beta);      
% Brerr=sqrt(gaussfilt_sph(Berrk_Pramsey.^2,gaz,gel,gaz,gel,sig_filt_beta));

%% VIS: fit params diagnostic
h=figure('Name','fit_param_diagnostic');
n_fit_param=size(Pramseyk_fpar,3);
for ii=1:n_fit_param
    ax=subplot(n_fit_param,1,ii);
    imagesc(rad2deg(az),rad2deg(el),Pramseyk_fpar(:,:,ii)');
    cbar=colorbar();
    title(cbar,Pramsey_cname{ii})
    ax.YDir='normal';
%     xlabel('$\theta$');
%     ylabel('$\phi$');
end

%% truncate sph-dist to within sensible elev-limits
% config
el_trunc=deg2rad(60);           % elevation angle threshold to truncate (rad)
% el_trunc=deg2rad(90);        
b_trunc=(abs(gel)>el_trunc);    % boolean indicator to truncate

% % get/store original
% % TODO - this needs work
% if ~exist('Bk_Pramsey_original','var')
%     % store
%     Bk_Pramsey_original=Bk_Pramsey;
%     Berrk_Pramsey_original=Berrk_Pramsey;
%     Br_original=Br;
%     Brerr_original=Brerr;
% else
%     % get
%     Bk_Pramsey=Bk_Pramsey_original;
%     Berrk_Pramsey=Berrk_Pramsey_original;
%     Br=Br_original;
%     Brerr=Brerr_original;
% end

% truncate
Bk_Pramsey(b_trunc)=NaN;
Berrk_Pramsey(b_trunc)=NaN;
% Br(b_trunc)=NaN;
% Brerr(b_trunc)=NaN;



%% Ramsey analysis: EQUATOR: P
% collate data to fit
yy=vertcat(P_eq{:});

Pramsey_fit_eq=fitnlm(1e6*tt,yy,Pramsey_mdl,...
    Pramsey_par0,'CoefficientNames',Pramsey_cname);

% get fit params
Pramsey_fpar_eq=Pramsey_fit_eq.Coefficients.Estimate;
Pramsey_fparerr_eq=Pramsey_fit_eq.Coefficients.SE;
    
% get freq
om_Pramsey_eq=Pramsey_fpar_eq(idx_om_Pramsey);
omerr_Pramsey_eq=Pramsey_fparerr_eq(idx_om_Pramsey);

% magnetic field
B_Pramsey_eq=1e6*om_Pramsey_eq/(2*pi*C_gymag);
Berr_Pramsey_eq=1e6*omerr_Pramsey_eq/(2*pi*C_gymag);


%% regions to display Ramsey signal
% max, min and (0,0)
loc_disp=[];    % row-array of (Iaz,Iel)-location to display
azel_disp=[];   % (az,el)

%%% theta,phi=(0,0)
% loc_disp=cat(1,loc_disp,[iaz_0,iel_0]);
% azel_disp=cat(1,azel_disp,[az(iaz_0),el(iel_0)]);
loc_disp(end+1,:)=[iaz_0,iel_0];
azel_disp(end+1,:)=[az(iaz_0),el(iel_0)];


%%% max -----------------------------------------------------------
% % equator
% [Bk_eq_max,iaz_Bmax]=max(Bk_eq);
% loc_disp=cat(1,loc_disp,[iaz_Bmax,iel_0]);
% azel_disp=cat(1,azel_disp,[az(iaz_Bmax),el(iel_0)]);

% anywhere
Bk_max=max(Bk_Pramsey(:));
[iaz_Bmax,iel_Bmax]=find(Bk_Pramsey==Bk_max);
loc_disp(end+1,:)=[iaz_Bmax,iel_Bmax];
azel_disp(end+1,:)=[az(iaz_Bmax),el(iel_Bmax)];

%%% min --------------------------------------------------------------
% % equator
% [Bk_eq_min,iaz_Bmin]=min(Bk_eq);
% loc_disp=cat(1,loc_disp,[iaz_Bmin,iel_0]);
% azel_disp=cat(1,azel_disp,[az(iaz_Bmin),el(iel_0)]);

% anywhere
Bk_min=min(Bk_Pramsey(:));
[iaz_Bmin,iel_Bmin]=find(Bk_Pramsey==Bk_min);
% loc_disp=cat(1,loc_disp,[iaz_Bmin,iel_Bmin]);
% azel_disp=cat(1,azel_disp,[az(iaz_Bmin),el(iel_Bmin)]);
loc_disp(end+1,:)=[iaz_Bmin,iel_Bmin];
azel_disp(end+1,:)=[az(iaz_Bmin),el(iel_Bmin)];

% summary
n_loc_disp=size(loc_disp,1);        % num of locations to display


%% VIS (publication): Mode-resolved Ramsey fringe: P
% markers to distinguish region
c_loc=zeros(n_loc_disp,3);
cl_loc=[0,1,1]'*[1,1,1];


figname='Pramsey_fringe_mode';
h=figure('Name',figname,'Units',config_fig.units,'Position',[0,0,8.6,4],'Renderer',config_fig.rend);
hold on;

pleg=[];
for ii=1:n_loc_disp
    iazel=loc_disp(ii,:);
    tazel=azel_disp(ii,:);
    
    tp=squeeze(P_k_avg(:,iazel(1),iazel(2)));
    tperr=squeeze(P_k_se(:,iazel(1),iazel(2)));
    
    dY_data=(ii-2)*0.5;
    
    pexp=ploterr(1e6*tau,dY_data + tp,[],tperr,config_fig.mark_typ{ii},'hhxy',0);
    set(pexp(1),'MarkerEdgeColor',c_loc(ii,:),...
        'MarkerFaceColor',cl_loc(ii,:),...
        'MarkerSize',config_fig.mark_siz,...
        'DisplayName',num2str(rad2deg(tazel)));
    set(pexp(2),'Color',c_loc(ii,:));
    pleg(ii)=pexp(1);
    
    yy=feval(Pramseyk_fit{iazel(1),iazel(2)},tt_fit);
    pfit=plot(tt_fit,dY_data + yy,...
        'LineWidth',config_fig.line_wid,'LineStyle',config_fig.line_sty{ii},'Color',c_loc(ii,:));
    uistack(pfit,'bottom');
end
ax=gca;
box on;

axis tight;
ax_xlim=ax.XLim;
ax_ylim=ax.YLim;

% ax.XLim=[1.9,3.6];
% ax.YLim=1.1*[-1,1];
xlim(ax_xlim+0.05*diff(ax_xlim)*[-1,1]);
ylim(ax_ylim+0.05*diff(ax_ylim)*[-1,1]);

set(ax,'Layer','Top');
ax.FontSize=config_fig.ax_fontsize;
ax.LineWidth=config_fig.ax_lwid;

xlabel('pulse delay $\tau~(\mu s)$');
ylabel('polarisation');

% lgd=legend(pleg);
% title(lgd,'$(\theta,~\phi)$ (deg)');

%% VIS (publication): Magnetometry (detected)
figname='halo_magnetometry_ff';
h=figure('Name',figname,'Units',config_fig.units,'Position',[0,0,8.6,4.5],'Renderer',config_fig.rend);

tp=plotFlatMapWrappedRad(gaz,gel,Bk_Pramsey,'rect','texturemap');

% for rect projection, IMAGESC --> controllable x,y axis?

% label ROI
hold on;
for ii=1:n_loc_disp
    tazel=rad2deg(azel_disp(ii,:));
    tp=plot(tazel(1),tazel(2),...
        'MarkerEdgeColor',c_loc(ii,:),'MarkerFaceColor',cl_loc(ii,:),...
        'Marker',config_fig.mark_typ{ii},'MarkerSize',config_fig.mark_siz);
end

% annotation
ax=gca;
set(ax,'Layer','Top');
box on;
% grid on;
ax.FontSize=config_fig.ax_fontsize;
ax.LineWidth=config_fig.ax_lwid;

% axis tight;
xlim([-180,180]);
ylim([-90,90]);
xticks(-180:90:180);
yticks(-90:45:90);

xlabel('$\theta$ (deg)');
ylabel('$\phi$ (deg)');

colormap('viridis');

%%% colorbar
% cbar=colorbar('west');
cbar=colorbar('eastoutside');
cbar.TickLabelInterpreter='latex';
% cbar.Label.Interpreter='latex';
% cbar.Label.String='$\mathrm{B}$ (G)';
% cbar.Label.FontSize=config_fig.ax_fontsize;
title(cbar,'$\mathrm{B}$ (G)','Interpreter','latex');       % title on top
cbar.FontSize=config_fig.ax_fontsize;
% change colorbar position
ax_pos = plotboxpos(ax);
pos_ax=get(gca,'Position');
pos_cbar=get(cbar,'Position');
pos_cbar(3)=0.025;
pos_cbar(1)=ax_pos(1)+ax_pos(3)+1*pos_cbar(3);
set(cbar,'Position',pos_cbar);
set(gca,'Position',pos_ax);     % return axis to original
% colorbar limits
Btomo_cbar_lim = cbar.Limits;
set(cbar,'Ticks',Btomo_cbar_lim);       % ticks ONLY at colorbar lims


% hatch-out truncated region --------------------------
% hack: two separate hatched regions
patch_xdata{1}=[ax.XLim(1), ax.XLim(2), ax.XLim(2), ax.XLim(1)];
patch_ydata{1}=[ax.YLim(1), ax.YLim(1), -rad2deg(el_trunc), -rad2deg(el_trunc)];

patch_xdata{2}=[ax.XLim(1), ax.XLim(2), ax.XLim(2), ax.XLim(1)];
patch_ydata{2}=[rad2deg(el_trunc), rad2deg(el_trunc), ax.YLim(2), ax.YLim(2)];

for ii=1:2
    hpatch_trunc(ii)=patch('XData',patch_xdata{ii},...
        'YData',patch_ydata{ii},...
        'FaceColor','none','EdgeColor','none');
    H_trunc(ii)=hatchfill2(hpatch_trunc(ii),'single','HatchAngle',45,'HatchDensity',20,...
        'HatchColor','k','HatchLineWidth',config_fig.line_wid);
    uistack(H_trunc(ii),'bottom');
end

% hpatch_trunc(2)=patch('XData',180*[-1,1,1,-1],...
%     'YData',[45,45,90,90],...
%     'FaceColor','none','EdgeColor','none');
% H_trunc(2)=hatchfill2(hpatch_trunc(2),'single','HatchAngle',45,'HatchDensity',20,...
%     'HatchColor','k','HatchLineWidth',config_fig.line_wid);

% TODO: get vecrast working with code below
% hpatch_trunc=patch('XData',[ax.XLim(1),ax.XLim(2),ax.XLim(2),ax.XLim(1)],...
%     'YData',[ax.YLim(1),ax.YLim(1),ax.YLim(2),ax.YLim(2)],...
%     'FaceColor','none','EdgeColor','none');
% uistack(hpatch_trunc,'bottom');
% H_trunc=hatchfill2(hpatch_trunc,'single','HatchAngle',45,'HatchDensity',20,...
%     'HatchColor','k','HatchLineWidth',config_fig.line_wid);


%---------------------------------------------
% NOTE: save with vecrast
% e.g.
% vecrast(h,strcat('B_dist_',getdatetimestr),600,'bottom','pdf')
%---------------------------------------------

%% VIS: histogram of B
figname='B_histogram_spatial';
h=figure('Name',figname,'Units',config_fig.units,'Position',[0,0,8.6,4],'Renderer',config_fig.rend);

hold on;

% naive one without lat-lon grid weights
X = Bk_Pramsey(~isnan(Bk_Pramsey));       % rid of NaNs
Bhist = histogram(X,'Normalization','pdf');
Bhist.DisplayStyle='stairs';
Bhist.EdgeColor='k';
Bhist.FaceColor='none';
Bhist.DisplayName='data';

% fit Gaussian
Xmean = mean(X(:));
Xstd = std(X(:));
xlim0=xlim;
xx = linspace(xlim0(1),xlim0(2),1e3);
yy = gaus_pdf(xx,Xmean,Xstd);
p_fit = plot(xx,yy,'r-','LineWidth',1);
p_fit.DisplayName='$\mathcal{N}(\bar\mu,\bar\sigma^2)$';

% set xlims to colorbar of tomography
xlim(Btomo_cbar_lim);

% annotate
ax=gca;
set(ax,'Layer','Top');
box on;
ax.FontSize=config_fig.ax_fontsize;
ax.LineWidth=config_fig.ax_lwid;
lgd=legend([Bhist,p_fit]);
lgd.Location='best' ;

xlabel('$\mathrm{B}$ (G)');
ylabel('spatial pdf');


%% VIS (publication): Magnetic tomography: equatorial: P 
figname='B_equatorial_tomography_P';
h=figure('Name',figname,'Units',config_fig.units,'Position',[0,0,8.6,4],'Renderer',config_fig.rend);
hold on;

% %%% EQ-integrated 
% H_eq_int=shadedErrorBar([-180,180],B_Pramsey_eq*[1,1],Berr_Pramsey_eq*[1,1],'k');
% H_eq_int.mainLine.LineWidth=1;
% H_eq_int.patch.FaceAlpha=config_fig.pf_alpha;
% H_eq_int.mainLine.DisplayName='$\theta$ int';
% H_eq_int.edge(1).Visible='off';
% H_eq_int.edge(2).Visible='off';


% %%% entire halo integrated (with BEC)
% % TODO - maybe elminate BECs (saturation) but this could be a point of
% % interesting discussion
% H_r_int=shadedErrorBar([-180,180],B_Pramsey_halo*[1,1],Berr_Pramsey_halo*[1,1],'k');
% H_r_int.mainLine.LineWidth=1;
% H_r_int.mainLine.LineStyle='--';
% H_r_int.mainLine.DisplayName='$\mathbf{r}$ int';
% H_r_int.patch.FaceAlpha=1;
% H_r_int.edge(1).Visible='off';
% H_r_int.edge(2).Visible='off';


%%% halo integrated (no BEC) 
H_halo_int=shadedErrorBar([-180,180],B_Pramsey_halo*[1,1],Berr_Pramsey_halo*[1,1],'k');
H_halo_int.mainLine.LineWidth=1;
H_halo_int.mainLine.LineStyle='--';
H_halo_int.mainLine.DisplayName='$\mathbf{r}$ int';
H_halo_int.mainLine.Visible='off';
H_halo_int.patch.FaceColor=0.8*ones(1,3);
H_halo_int.patch.FaceAlpha=1;
H_halo_int.edge(1).Visible='off';
H_halo_int.edge(2).Visible='off';


%%% spatial resolved predictions
%%%%% Far-Field
% %PLOTERR
% tp=ploterr(rad2deg(az),Bk_eq,[],Bkerr_eq,'o','hhxy',0);
% set(tp(1),'Marker','.','MarkerSize',4.5,...
%     'MarkerFaceColor',config_fig.coll_theme(2,:),'MarkerEdgeColor',config_fig.col_theme(2,:),...
%     'DisplayName','');
% set(tp(2),'Color',config_fig.col_theme(2,:));

% SHADED ERR BAR
H_res_ff=shadedErrorBar(rad2deg(az),Bk_eq,Bkerr_eq,'r');
H_res_ff.mainLine.Color=config_fig.col_theme(2,:);
H_res_ff.mainLine.LineStyle='-';
H_res_ff.mainLine.LineWidth=1;
H_res_ff.mainLine.DisplayName='$\mathbf{r}$ resolved $\infty$';
H_res_ff.patch.FaceColor=config_fig.coll_theme(2,:);  
H_res_ff.patch.FaceAlpha=config_fig.pf_alpha;
% H_res_ff.edge(1).Visible='off';
% H_res_ff.edge(2).Visible='off';
H_res_ff.edge(1).Color=config_fig.coll_theme(2,:);
H_res_ff.edge(2).Color=config_fig.coll_theme(2,:);


% %%% ROI
% for ii=1:n_loc_disp
%     clearvars tp;
%     
%     iazel=loc_disp(ii,:);
%     taz_disp=azel_disp(ii,1);
%     tp=ploterr(rad2deg(taz_disp),Bk_eq(iazel(1)),[],Bkerr_eq(iazel(1)),...
%         config_fig.mark_typ{ii},'hhxy',1);
%     
%     set(tp(1),'MarkerEdgeColor',c_loc(ii,:),...
%         'MarkerFaceColor',cl_loc(ii,:),...
%         'MarkerSize',config_fig.mark_siz,...
%         'DisplayName',num2str(rad2deg(taz_disp)));
%     set(tp(2),'Color',c_loc(ii,:));
%     
% end


%%% annotation
box on;
ax=gca;
set(ax,'Layer','Top');
ax.FontSize=config_fig.ax_fontsize;
ax.LineWidth=config_fig.ax_lwid;

xlabel('$\theta$ (deg)');
ylabel('$\mathrm{B}$ (G)');

ax.XTick=-180:90:180;
axis_snug(ax,[0,0.1]);

% % legend
% lgd=legend([H_res_r.mainLine,H_res_ff.mainLine,H_r_int.mainLine]);
% lgd.Location='eastoutside';


% %%% Zoomed (fluctuation) - INSET -------------------------------------
% ax2=axes('Position',[.7 .5 .2 .4]);      % 'YAxisLocation','right','XAxisLocation','top'
% set(ax2,'Layer','top');
% box on;
% hold on; 
% 
% % Ramsey
% %%% halo integrated (no BEC) 
% H_halo_int=shadedErrorBar([-180,180],B_Pramsey_halo*[1,1],Berr_Pramsey_halo*[1,1],'k');
% H_halo_int.mainLine.LineWidth=1;
% H_halo_int.mainLine.LineStyle='--';
% H_halo_int.mainLine.DisplayName='$\mathbf{r}$ int';
% H_halo_int.mainLine.Visible='off';
% H_halo_int.patch.FaceColor=0.8*ones(1,3);
% H_halo_int.patch.FaceAlpha=1;
% H_halo_int.edge(1).Visible='off';
% H_halo_int.edge(2).Visible='off';
% 
% %%% spatial resolved predictions
% H_res_ff=shadedErrorBar(rad2deg(az),Bk_eq,Bkerr_eq,'r');
% H_res_ff.mainLine.Color=config_fig.col_theme(2,:);
% H_res_ff.mainLine.LineStyle='-';
% H_res_ff.mainLine.LineWidth=1;
% H_res_ff.mainLine.DisplayName='$\mathbf{r}$ resolved $\infty$';
% H_res_ff.patch.FaceColor=config_fig.coll_theme(2,:);  
% H_res_ff.patch.FaceAlpha=config_fig.pf_alpha;
% H_res_ff.edge(1).Color=config_fig.coll_theme(2,:);
% H_res_ff.edge(2).Color=config_fig.coll_theme(2,:);
% 
% % annotation
% % zoom to theta=0
% xlim(10*[-1,1]);
% ylim([0.528,0.545]);
% ax2.FontSize=ax.FontSize-1;
% ax2.TickLength=2*ax2.TickLength;    % increase ticklength


%% BB B-difference around equator
dth=mean(diff(az));     % incremental diff angle scanned around equator (rad)
dI_pi=round(pi/dth);    % # increments to shift for azimuthal BB

az_bb=circshift(az,dI_pi);      % azim-vectors pi-shifted

% check BB
dbb_max=max(abs(abs(az_bb-az)-pi));     % max deviation from BB (rad)
if dbb_max>1e-2
    warning('Some regions are not back-to-back!');
end

% get B-params at opposite regions
Bk_eq_bb=circshift(Bk_eq,dI_pi);
Bkerr_eq_bb=circshift(Bkerr_eq,dI_pi);

% BB-difference
% dB_eq_bb=abs(Bk_eq_bb - Bk_eq);
dB_eq_bb=Bk_eq - Bk_eq_bb;
dBerr_eq_bb=sqrt(Bkerr_eq_bb.^2 + Bkerr_eq.^2);     % uncorrelated noise

%%% vis
figname='dBbb_ramsey_ff';
h=figure('Name',figname,'Units',config_fig.units,'Position',[0,0,8.6,6],'Renderer',config_fig.rend);

% B-field around equator @ theta and theta+pi
subplot(2,1,1);
hold on;
% plot(rad2deg(az),Bk_eq);            % B-field (G)
% plot(rad2deg(az),Bk_eq_bb);         % B-field at BB region (G)
clearvars tp;
tp(1)=shadedErrorBar(rad2deg(az),Bk_eq,Bkerr_eq,'r',1);        % B-field (G)
tp(2)=shadedErrorBar(rad2deg(az),Bk_eq_bb,Bkerr_eq_bb,'b',1);  % B-field at BB region (G)

ax=gca;
set(ax,'Layer','Top');
xlim([0,180]);      % periodic
xticks(0:45:180);
xlabel('$\theta$ (deg)');
ylabel('$B$ (G)');

ax.FontSize=9;
% ax.LineWidth=1;

% lgd=legend([tp(1).mainLine,tp(2).mainLine],{'$\theta$','$\theta+\pi$'});

% B-field difference @ BB regions around equator
subplot(2,1,2);
hold on;

tp=shadedErrorBar(rad2deg(az),1e3*dB_eq_bb,1e3*dBerr_eq_bb,'k');
% tp.mainLine.Color=config_fig.col_theme(2,:);
% tp.mainLine.LineWidth=1;
% tp.patch.FaceColor=config_fig.coll_theme(2,:);
% tp.patch.FaceAlpha=config_fig.pf_alpha;
% tp.edge(1).Visible='off';
% tp.edge(2).Visible='off';

ax=gca;
set(ax,'Layer','Top');
xlim([0,180]);      % periodic
ylim_0=ax.YLim;
% ylim([0,ylim_0(2)]);
xticks(0:45:180);
xlabel('$\theta$ (deg)');
ylabel('$\Delta B_{\mathrm{BB}}$ (mG)');

ax.FontSize=9;
% ax.LineWidth=1;


%% dBdx (equator)
% mean estimate of dBdx
dBdx_eq=dB_eq_bb/d_sep;     % in G/m

% standard error
%       in progress
% error in separation dist
dx_frerr=sqrt(2)*sig_psf_r/d_sep;       % frac-err of diff X
dB_frerr=dBerr_eq_bb./dB_eq_bb;         % frac-err of diff B
dBdx_eq_frerr=sqrt(dx_frerr^2+dB_frerr.^2); % frac-err of dBdx
dBdx_eq_se=abs(dBdx_eq_frerr.*dBdx_eq);

% plot --------------------------------------------
figname='dBdx_ramsey_ff';
h=figure('Name',figname,'Units',config_fig.units,'Position',[0,0,8.6,6],'Renderer',config_fig.rend);
tp=shadedErrorBar(rad2deg(az),dBdx_eq,dBdx_eq_se,'k');

ax=gca;
set(ax,'Layer','Top');

xlim([0,180]);      % periodic
ylim_Original=ax.YLim;
% ylim([0,ylim_Original(2)]);

xticks(0:45:180);

xlabel('$\theta$ (deg)');
ylabel('$d\mathrm{B}/dx$ (G/m)');


%% save outputs
vars_to_save={'az','el','gaz','gel','alpha',...
    'configs','config_fig',...
    'iaz_0','iel_0','n_loc_disp','azel_disp',...
    'tau','P_k_avg','P_k_se',...
    'Pramseyk_fit','P_k',...
    'tt_fit',...
    'B_Pramsey_halo','Berr_Pramsey_halo',...
    'Bk_Pramsey','Berrk_Pramsey',...
    'Bk_eq','Bkerr_eq',...
    'dBdx_eq','dBdx_eq_se'};  

% check exist
for ii=1:numel(vars_to_save)
    tvarname=vars_to_save{ii};
    if ~exist(tvarname,'var')
        warning('variable %s does not exist.',tvarname);
    end
end

save(['out_',getdatetimestr,'.mat'],vars_to_save{:});

%% end of script
toc