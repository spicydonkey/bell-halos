% Characterising global rotation around Y-axis
%
%
% DKS
% 2018-05-29

%% CONFIGS
%%% David PC
% config_name='C:\Users\David\Documents\MATLAB\bell-halos\analysis\exp3_yrot_char\src\config_1.m';
% config_name='C:\Users\David\Documents\MATLAB\bell-halos\analysis\exp3_yrot_char\src\config_2.m';
% config_name='C:\Users\David\Documents\MATLAB\bell-halos\analysis\exp3_yrot_char\src\config_s1.m';

%%% HE BEC PC
config_name='C:\Users\HE BEC\Documents\MATLAB\bell-halos\analysis\exp3_yrot_char\src\config_s1_halocapture.m';
% config_name='C:\Users\HE BEC\Documents\MATLAB\bell-halos\analysis\exp3_yrot_char\src\config_s1.m';
% config_name='C:\Users\HE BEC\Documents\MATLAB\bell-halos\analysis\exp3_yrot_char\src\config_1.m';
% config_name='C:\Users\HE BEC\Documents\MATLAB\bell-halos\analysis\exp3_yrot_char\src\config_2.m';

% vis
f_units='normalized';
f_pos=[0.2,0.2,0.2,0.3];
f_ren='painters';

[c,cl,cd]=palette(3);
% c_gray=0.6*ones(1,3);
% line_sty={'-','--',':'};
mark_typ={'o','^','d','s'};
% str_ss={'$\vert\!\uparrow\rangle$','$\vert\!\downarrow\rangle$'};
% str_ss={'$m_J = 1$','$m_J = 0$'};
% str_ss={'$m_J = 1$','$m_J = 0$','$m_J = -1$'};
% str_ss={'$\uparrow\uparrow$','$\downarrow\downarrow$','$\uparrow\downarrow$'};
mark_siz=7;
line_wid=1.5;
fontsize=12;
ax_lwidth=1.2;

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
    % TODO
    %   do I need to set some things to default? or re-code analysis?
    
    % defaults autoscan-related vars
    nparam=1;       % 1- since we don't search any params
    id_in_param={configs.load.id};    % all IDS to load
    
    par_T=NaN;      % default param to NaN
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



%% distinguish mJ and capture halo
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
%%% unit spherise
k_halo=zxy0_filt;       % initialise atoms in k-space (i.e. atoms lie on unit-sphere)

for ii=1:n_mf
    k_halo(:,ii)=cellfun(@(x) x/r_halo_avg(ii),k_halo(:,ii),'UniformOutput',false);
end

%%% ellipsoid fit to sphere
for ii=1:n_mf
    k_halo(:,ii)=map2usph(k_halo(:,ii));
end
    
% DEBUG
scatter_halo(k_halo);



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



%% categorise data by exp-params
k_par=cell(1,nparam);
if nparam>1
    for ii=1:nparam
        k_par{ii}=k_halo_filt(b_paramset(:,ii),:);      % get all halo data
        %from this param-set and store
    end
else
    k_par{1}=k_halo_filt;
end

% DEBUG
figure;
for ii=1:nparam
    subplot(1,nparam,ii);
    plot_zxy(k_par{ii});
    
    axis equal;
    xlabel('kx');
    ylabel('ky');
    zlabel('kz');
    view([0,0]);
end


%%% END OF PREPROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% clean workspace
clear txy zxy zxy0 zxy0_filt tzxy tzxy_mf tp_bec tp_bec0 tp_halo tr_bec0 tr_lim;
clear h_zxy*;       % clear figs

%% ANALYSIS

k_par_orig=k_par;       % save original

%% HACK: optimise halo centering
k_par=k_par_orig;       % to orig

% displace
for ii=1:nparam
    tk=k_par{ii};
    
    for jj=1:n_mf
        tk(:,jj)=boost_zxy(tk(:,jj),configs.post.Dk{jj});
    end
    
    k_par{ii}=tk;
end

%% non-dimensionalzed time
T_larmor=0.7e-6;        % Larmor precession period [s]
tau_rot=par_T/T_larmor;


%% Number counting - momentum-unresolved
%   NOTE: Nsc Poissonian and stat analysis needs to be careful

latlon_halo_ed={[-pi,pi],[-pi/2,pi/2]};     % momentum integrated
nsc_halo=cell(1,nparam);
for ii=1:nparam
    nsc_halo{ii}=cell(1,n_mf);
    temp_nsc_halo=cellfun(@(k) histlatlon_halo(k,latlon_halo_ed),k_par{ii},'UniformOutput',false);
    
    % tidy form by collapsing independent shots into dim-3
    for jj=1:n_mf
        nsc_halo{ii}{jj}=cat(3,temp_nsc_halo{:,jj});
    end
end

%%% mJ population fraction per mode/shot
%   NOTE: this mode is the whole halo
% sum scattered num. in each (zone, mJ) / shot
nsctot_halo=cellfun(@(n) sum(cat(4,n{:}),4),nsc_halo,'UniformOutput',false);    % total number per mode/shot

% population fraction
P_mJ_halo=cell(1,nparam);
for ii=1:nparam
    temp_nsc_tot=nsctot_halo{ii};       % total scattered num array for this subset
    for jj=1:n_mf
        P_mJ_halo{ii}{jj}=nsc_halo{ii}{jj}./temp_nsc_tot;   % eval pop fraction
    end
end

% statistics
temp_P_mJ_halo_avg=cell(1,nparam);
temp_P_mJ_halo_std=cell(1,nparam);
for ii=1:nparam
    temp_P_mJ_halo_avg{ii}=cellfun(@(p) mean(p,3),P_mJ_halo{ii},'UniformOutput',false);
    temp_P_mJ_halo_std{ii}=cellfun(@(p) std(p,[],3),P_mJ_halo{ii},'UniformOutput',false);
end
% tidy form
P_mJ_halo_avg=cell(1,n_mf);
P_mJ_halo_std=cell(1,n_mf);
temp_P_mJ_halo_avg=cat(1,temp_P_mJ_halo_avg{:});      % collapse param# to cell-array row#
temp_P_mJ_halo_std=cat(1,temp_P_mJ_halo_std{:});

% collapse param# as dim-3 in 2D histogram for each mJ
for ii=1:n_mf
    P_mJ_halo_avg{ii}=cat(3,temp_P_mJ_halo_avg{:,ii});
    P_mJ_halo_std{ii}=cat(3,temp_P_mJ_halo_std{:,ii});
end

% 1x1 grid to squeeze into array
for ii=1:n_mf
    P_mJ_halo_avg{ii}=squeeze(P_mJ_halo_avg{ii});
    P_mJ_halo_std{ii}=squeeze(P_mJ_halo_std{ii});
end
P_mJ_halo_avg=cat(2,P_mJ_halo_avg{:});
P_mJ_halo_std=cat(2,P_mJ_halo_std{:});


%% Fit Rabi oscillation (momentum-modes unresolved)
%   General near-resonant Rabi oscillation model for a TLS: text
%   representation

tt=linspace(0,max(par_T),1e3);

% (A) simple phenomenological (2-param): 
%   * amplitude (pkpk) [1]
%   * Rabi frequency [rad/s]
rabi_mdl{1}='p~1-0.5*amp*(1-cos(om*x1))';    % for the initially unpopulated state mJ=1
rabi_mdl{2}='p~0.5*amp*(1-cos(om*x1))';    % for the initially unpopulated state mJ=0


fit_rabi_halo=cell(2,1);
fit_rabi_halo_params=cell(2,1);
pp=cell(2,1);
P_errfrac=cell(2,1);

% idx_mJ=2;    % mJ=0
for idx_mJ=1:2
    tp=P_mJ_halo_avg(:,idx_mJ);
    
    % estimate params
    tRabiAmp=max(tp)-min(tp);
    tRabiOmega=2*pi/20e-5;      % we keep this near const
    tparam0=[tRabiAmp,tRabiOmega];
    
    %%% fit model
    fopts_halo=statset('Display','off');
    fit_rabi_halo{idx_mJ}=fitnlm(par_T,tp,rabi_mdl{idx_mJ},tparam0,'CoefficientNames',{'amp','om'},...
        'Options',fopts_halo);
    fit_rabi_halo_params{idx_mJ}=fit_rabi_halo{idx_mJ}.Coefficients.Estimate;   % get fitted params
    
    % evaluate fitted model
    pp{idx_mJ}=feval(fit_rabi_halo{idx_mJ},tt);
    
    %%% Uncertainty model
    %   Constant error fraction
    P_errfrac{idx_mJ}=harmmean(P_mJ_halo_std(:,idx_mJ)./P_mJ_halo_avg(:,idx_mJ));      % gives a reasonal fit to ~pi/2
end


%% MOMENTUM-RESOLVED: mJ-oscillation
%   
%   1. define momentum "zones" (multiple momentum modes)
%       e.g.    (A) lat-lon zones
%               (B) conical section of halo: orientation-vector and inclusion half-angle
%   2. for each halo/exp#/expparam get #atoms in each zone
%   3. compare for each momentum zone/expparam do statistics for n_mJ/n_tot
%   4. zone resolved Rabi oscillation
%   5. fit oscillation characteristic per momentum zone --> e.g. resonant freq + detuning model
%
%   * uniformity study
%   * structure around halo?
%   * extrapolate to momentum mode? dependency on zone size?
%

%%% define momentum zones
% (A) lat-lon grid
nzone_th=4;     % num. zones to equipartition azimuthal (-pi,pi]
nzone_phi=2;    % for elevation [-pi/2,pi/2]
momzone_th=linspace(-pi,pi,nzone_th+1);
% momzone_phi=linspace(-pi/2,pi/2,nzone_phi+1);
momzone_phi=linspace(-0.8,0.8,nzone_phi+1);     % elev limits to BEC caps

% TODO
% (B) conical zone (either overlaps or misses)


%%% histogram atoms into zones
% (A) lat-lon grid: 2D histogram
latlon_zone_ed={momzone_th,momzone_phi};
nsc_zone=cell(1,nparam);
for ii=1:nparam
    nsc_zone{ii}=cell(1,n_mf);
    temp_nsc_zone=cellfun(@(k) histlatlon_halo(k,latlon_zone_ed),k_par{ii},'UniformOutput',false);

    % tidy form by collapsing independent shots into dim-3
    for jj=1:n_mf
        nsc_zone{ii}{jj}=cat(3,temp_nsc_zone{:,jj});
    end
end


%%% mJ population fraction per mode/shot
% sum scattered num. in each (zone, mJ) / shot
nsctot_zone=cellfun(@(n) sum(cat(4,n{:}),4),nsc_zone,'UniformOutput',false);    % total number per mode/shot

% population fraction
P_mJ_zone=cell(1,nparam);
for ii=1:nparam
    temp_nsc_tot=nsctot_zone{ii};       % total scattered num array for this subset
    for jj=1:n_mf
        P_mJ_zone{ii}{jj}=nsc_zone{ii}{jj}./temp_nsc_tot;   % eval pop fraction
    end
end

% statistics
temp_P_mJ_zone_avg=cell(1,nparam);
temp_P_mJ_zone_std=cell(1,nparam);
for ii=1:nparam
    temp_P_mJ_zone_avg{ii}=cellfun(@(p) mean(p,3),P_mJ_zone{ii},'UniformOutput',false);
    temp_P_mJ_zone_std{ii}=cellfun(@(p) std(p,[],3),P_mJ_zone{ii},'UniformOutput',false);
end

% tidy form
P_mJ_zone_avg=cell(1,n_mf);
P_mJ_zone_std=cell(1,n_mf);

temp_P_mJ_zone_avg=cat(1,temp_P_mJ_zone_avg{:});      % collapse param# to cell-array row#
temp_P_mJ_zone_std=cat(1,temp_P_mJ_zone_std{:});

% collapse param# as dim-3 in 2D histogram for each mJ
for ii=1:n_mf
    P_mJ_zone_avg{ii}=cat(3,temp_P_mJ_zone_avg{:,ii});
    P_mJ_zone_std{ii}=cat(3,temp_P_mJ_zone_std{:,ii});
end



%% Rabi oscillation (momentum-zone resolved fit)
idx_mJ=2;   % ONLY mJ=0

% momzone_amp=NaN(nzone_th,nzone_phi);
% momzone_om=NaN(nzone_th,nzone_phi);

P_momzone_0=P_mJ_zone_avg{idx_mJ};  % get pop fracs for all (expparams,zones) for mJ=0

fit_rabi_zone=cell(nzone_th,nzone_phi);
fit_rabi_zone_param=NaN(nzone_th,nzone_phi,2);      % stores amp/omega in dim3

for ii=1:nzone_th*nzone_phi
    [mm,nn]=ind2sub([nzone_th,nzone_phi],ii);  % get this zone
    tp=squeeze(P_momzone_0(mm,nn,:));    % pop fraction profile for this zone
    
    % check if this zone has any atoms
    if sum(isnan(tp))>0.33*length(tp)     % all NaN (a little too strict)
        continue
    end
        
    % estimate params
    tRabiAmp=max(tp)-min(tp);
    tRabiOmega=2*pi/(20e-5);        % we keep this near const
    
    tparam0=[tRabiAmp,tRabiOmega];
    
    % fit model
    fopts_raw=statset('Display','off');
    fit_rabi_zone{mm,nn}=fitnlm(par_T,tp,rabi_mdl{idx_mJ},tparam0,'CoefficientNames',{'amp','om'},...
        'Options',fopts_raw);
    fit_rabi_zone_param(mm,nn,:)=fit_rabi_zone{mm,nn}.Coefficients.Estimate;
end
% get statistics around halo
rabi_amp_zone=fit_rabi_zone_param(:,:,1);
errfrac_amp=std(rabi_amp_zone(:))/mean(rabi_amp_zone(:));



%% DATA VISUALIZATION
%% Momentum mode unresolved
%%% plot
h_rabi_halo=figure('Name','rabi_halo','Units',f_units,'Position',[0.2,0.2,0.2,0.25],'Renderer',f_ren);
hold on;

%%% model
% % Fitted model to avg (a composite graphics object) - plotted first to be in bottom
% lineProps.col={cl(idx_mJ,:)};
% tfitted=mseb(1e6*tt,pp,P_errfrac.*pp,lineProps,1);

% momentum zone variability
for ii=1:2
    pp=feval(fit_rabi_halo{ii},tt);
    lineProps.col={cl(ii,:)};
    
    if ii==1
        p_err=errfrac_amp.*(1-pp);
    elseif ii==2
        p_err=errfrac_amp.*pp;
    end
    tfitted=mseb(1e6*tt,pp,p_err,lineProps,1);
end

%%% data
h=NaN(n_mf,1);
for ii=1:n_mf
    th=ploterr(1e6*par_T,P_mJ_halo_avg(:,ii),[],P_mJ_halo_std(:,ii),'o','hhxy',0);
%     th=ploterr(tau_rot,P_mJ_halo_avg(:,ii),[],P_mJ_halo_std(:,ii),'o','hhxy',0);
    set(th(1),'color',c(ii,:),'Marker',mark_typ{ii},'LineWidth',line_wid,...
        'MarkerSize',mark_siz,'MarkerFaceColor',cl(ii,:),...
        'DisplayName',num2str(configs.mf(ii).mf));
    set(th(2),'color',c(ii,:),'LineWidth',line_wid);
    
    h(ii)=th(1);
end

% annotations
box on;
ax=gca;
ax.LineWidth=ax_lwidth;
ax.FontSize=fontsize;
set(ax,'Layer','Top');     % graphics axes should be always on top
axis tight;

ylim([0,1]);
% ax.YTick=0:0.2:1;

xlabel('Pulse duration [$\mu$s]');
% xlabel('Pulse duration, $\tau$');
% ylabel('$P$');
ylabel('Population fraction $P$');

lgd=legend(h,'Location','East');
% title(lgd,'$m_J$');
set(lgd,'FontSize',fontsize-1);

h=gcf;
h.Renderer=f_ren;

% %% raw data (non-averaged)
% %%% collate all the data
% % pulse duration
% % T_collate=cellfun(@(x) ones(size(k,1),1),k_par,'UniformOutput',false);
% T_collate=cell(nparam,1);
% for ii=1:nparam
%     T_collate{ii}=par_T(ii)*ones(size(k_par{ii},1),1);
% end
% T_collate=cat(1,T_collate{:});
% 
% % mJ=0 pop fraction
% P_collate=cellfun(@(p) squeeze(p{2}),P_mJ_halo,'UniformOutput',false);
% P_collate=cat(1,P_collate{:});
% 
% 
% %%% FIT - same as for avgs
% % estimate params
% tRabiAmp=max(P_collate)-min(P_collate);
% tRabiOmega=2*pi/20e-5;      % we keep this near const
% tparam0=[tRabiAmp,tRabiOmega];
% 
% %%% fit model
% fopts_raw=statset('Display','off');
% fit_rabi_halo_raw=fitnlm(T_collate,P_collate,rabi_mdl,tparam0,'CoefficientNames',{'amp','om'},...
%     'Options',fopts_raw);
% 
% % evaluate fitted model
% tt=linspace(0,max(par_T),1e3);
% pp=feval(fit_rabi_halo_raw,tt);
% 
% 
% %%% Data vis
% h=figure('Name','rabi_halo_raw');
% 
% % raw data
% plot(T_collate*1e6,P_collate,'r.');
% % fit
% hold on;
% plot(tt*1e6,pp,'r-');
% 
% % annotation
% ax=gca;
% axis tight;
% box on;
% 
% ylim([0,1]);
% ax.YTick=0:0.2:1;
% 
% xlabel('Pulse duration [$\mu$s]');
% % xlabel('Pulse duration, $\tau$');
% ylabel('$P$');
% 

%% Momentum resolved
%%% COLORSPACES for many lines
% 1. max. distinguishable
% [c0_zone,clight_zone,cdark_zone]=palette(nzone_th*nzone_phi);

% 2. smooth colormaps
% 2.1. parula (dark) - seems best
% clight_zone=parula(nzone_th*nzone_phi);
% [~,c0_zone]=colshades(clight_zone);

% 2.2. viridis (dark)
clight_zone=viridis(nzone_th*nzone_phi);
[~,c0_zone]=colshades(clight_zone);

% 2.3. plasma (dark) - looks like magma
% clight_zone=plasma(nzone_th*nzone_phi);
% [~,c0_zone]=colshades(clight_zone);

h_rabi_momzone=figure('Name','rabi_momzone','Units',f_units,'Position',f_pos,'Renderer',f_ren);
hold on;

h=[];
for ii=1:n_mf
    % iterate over lat-long grid 2D indices
    for jj=1:nzone_th*nzone_phi
        [mm,nn]=ind2sub([nzone_th,nzone_phi],jj);
        th=ploterr(1e6*par_T,squeeze(P_mJ_zone_avg{ii}(mm,nn,:)),[],...
            squeeze(P_mJ_zone_std{ii}(mm,nn,:)),'o','hhxy',0);
%         th=ploterr(tau_rot,squeeze(P_mJ_zone_avg{ii}(mm,nn,:)),[],...
%             squeeze(P_mJ_zone_std{ii}(mm,nn,:)),'o','hhxy',0);
        set(th(1),'color',c0_zone(jj,:),'Marker',mark_typ{ii},'LineWidth',line_wid,...
            'MarkerSize',mark_siz,'MarkerFaceColor',clight_zone(jj,:));      % 'DisplayName',num2str(configs.mf(ii).mf)
        set(th(2),'color',c0_zone(jj,:),'LineWidth',line_wid);
        
        h=cat(1,h,th(1));
    end
end

for ii=1:nzone_th*nzone_phi
    [mm,nn]=ind2sub([nzone_th,nzone_phi],ii);  % get this zone
    % evaluate fitted model in this zone
    pp=feval(fit_rabi_zone{mm,nn},tt);
    
    % vis
    figure(h_rabi_momzone);
    hold on;
    tfitted=plot(1e6*tt,pp,'-','Color',c0_zone(ii,:));
    uistack(tfitted,'bottom');
end
    
% annotations
box on;
ax=gca;
ax.LineWidth=ax_lwidth;
ax.FontSize=fontsize;
set(ax,'Layer','Top');     % graphics axes should be always on top

axis tight;
% ax.YTick=0:0.2:1;
ylim([0,1]);

xlabel('Pulse duration [$\mu$s]');
% xlabel('Pulse duration, $\tau$');
% ylabel('$P$');
ylabel('Population fraction $P$');

%% Exploded figure of halo by lat-lon zones
%   could be refactored as a function
%
%   1. collate k vecs
%   2. for each zone find included k
%   3. translate by zone's center vector q
%   4. shade/color-coded scatter plot
%

%%% collate momentum-vecs
k_collated=cat(1,k_par{:});         % ALL data
k_collated=cat(1,k_collated{:,2});    % ONLY mJ=0

ks_collated=zxy2sphpol(k_collated);     % sph-polar vecs
k_th=ks_collated(:,1);
k_phi=ks_collated(:,2);

%%% sort k-vecs into zones
k_momzone=cell(nzone_th,nzone_phi);
for ii=1:numel(k_momzone)
    [mm,nn]=ind2sub(size(k_momzone),ii);    % get zone

    % boolean checks for vec in this zone
    tth_in=(k_th>momzone_th(mm))&(k_th<momzone_th(mm+1));
    tphi_in=(k_phi>momzone_phi(nn))&(k_phi<momzone_phi(nn+1));
    tzone_in=tth_in&tphi_in;
    
    k_momzone{ii}=k_collated(tzone_in,:);
end

%%% explode zones radially
disp_expl=0.1;      % displacement of explosion

momzone_th_cent=momzone_th(1:end-1)+0.5*diff(momzone_th);
momzone_phi_cent=momzone_phi(1:end-1)+0.5*diff(momzone_phi);

k_momzone_exp=cell(nzone_th,nzone_phi);
for ii=1:numel(k_momzone_exp)
    [mm,nn]=ind2sub(size(k_momzone),ii);    % get zone
    
    % vector to zone-centre
    tdk_expl=sphpol2zxy([momzone_th_cent(mm),momzone_phi_cent(nn),disp_expl]);
    
    k_momzone_exp{ii}=k_momzone{ii}+tdk_expl;
end


%%% DATA VIS
plot_scat_size=1;

% colorcoding
% gray
plot_scat_col=gray(numel(k_momzone_exp)+2);
plot_scat_col=plot_scat_col(2:end-1,:);

% colormaps
% plot_scat_col=palette(numel(k_momzone_exp));


h_halo_exp=figure('Name','halo_momzone_exploded');
hold on;
for ii=1:numel(k_momzone_exp)
    [mm,nn]=ind2sub(size(k_momzone),ii);    % get zone
    
    scatter_zxy(k_momzone_exp{mm,nn},plot_scat_size,c0_zone(ii,:));
    % color to match the zonal Rabi plot
end

view(3);
axis equal;
axis off;
