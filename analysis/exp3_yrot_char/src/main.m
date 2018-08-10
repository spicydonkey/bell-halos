% Characterising global rotation around Y-axis
%
%
% DKS
% 2018-05-29

%%% David PC
% config_name='C:\Users\David\Documents\MATLAB\bell-halos\analysis\exp3_yrot_char\src\config_1.m';
% config_name='C:\Users\David\Documents\MATLAB\bell-halos\analysis\exp3_yrot_char\src\config_2.m';
% config_name='C:\Users\David\Documents\MATLAB\bell-halos\analysis\exp3_yrot_char\src\config_s1.m';

%%% HE BEC PC
config_name='C:\Users\HE BEC\Documents\MATLAB\bell-halos\analysis\exp3_yrot_char\src\config_s1.m';
% config_name='C:\Users\HE BEC\Documents\MATLAB\bell-halos\analysis\exp3_yrot_char\src\config_1.m';
% config_name='C:\Users\HE BEC\Documents\MATLAB\bell-halos\analysis\exp3_yrot_char\src\config_2.m';

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


%% Below analysis for mJ-population oscillation is averaged in momentum
%% number of counts captured in halo
%   NOTE: Nsc Poissonian and stat analysis needs to be careful

% number of counts detected in filtered halo region
n_sc_counts=cellfun(@(v) shotSize(v),k_par,'UniformOutput',false);

% statistics: avg and std for number in each halo
n_sc_counts_avg=cellfun(@(n) mean(n,1),n_sc_counts,'UniformOutput',false);
n_sc_counts_avg=cat(1,n_sc_counts_avg{:});

n_sc_counts_std=cellfun(@(n) std(n,[],1),n_sc_counts,'UniformOutput',false);
n_sc_counts_std=cat(1,n_sc_counts_std{:});
    

% output
fprintf('Summary: Counts in halo\n');
for ii=1:nparam
    fprintf('par %d',ii);
    for jj=1:n_mf
        fprintf('\t:\t%5.3g(%2.2g)',[n_sc_counts_avg(ii,jj),n_sc_counts_std(ii,jj)]);
    end
    fprintf('\n');
end


%% number transfer in the triplet
n_sc_tot=cellfun(@(n) sum(n,2),n_sc_counts,'UniformOutput',false);      % total sc counts det'd in each shot

% state population fraction and Rabi oscillation
P_rabi_shot=cellfun(@(n,N) n./N,n_sc_counts,n_sc_tot,'UniformOutput',false);
P_rabi_avg=cellfun(@(p) mean(p,1),P_rabi_shot,'UniformOutput',false);
P_rabi_avg=cat(1,P_rabi_avg{:});
P_rabi_std=cellfun(@(p) std(p,[],1),P_rabi_shot,'UniformOutput',false);
P_rabi_std=cat(1,P_rabi_std{:});


%% non-dimensionalzed time
T_larmor=0.7e-6;        % Larmor precession period [s]
tau_rot=par_T/T_larmor;


%% Fit Rabi oscillation (momentum-modes unresolved)
%   General near-resonant Rabi oscillation model for a TLS: text
%   representation

% (A) simple phenomenological (2-param): 
%   * amplitude (pkpk) [1]
%   * Rabi frequency [rad/s]
rabi_mdl='p~0.5*amp*(1-cos(om*x1))';    % for the initially unpopulated state mJ=0

halo_amp=NaN;
halo_om=NaN;

idx_mJ=2;    % mJ=0
tp=P_rabi_avg(:,idx_mJ);

% estimate params
tRabiAmp=max(tp)-min(tp);
tRabiOmega=2*pi/20e-5;      % we keep this near const
tparam0=[tRabiAmp,tRabiOmega];

%%% fit model
tfopts=statset('Display','off');
tfit_rabi=fitnlm(par_T,tp,rabi_mdl,tparam0,'CoefficientNames',{'amp','om'},...
    'Options',tfopts);
tparam_fit=tfit_rabi.Coefficients.Estimate;

% get fitted model params
halo_amp=tparam_fit(1);
halo_om=tparam_fit(2);

% evaluate fitted model
tt=linspace(0,max(par_T),1e3);
pp=feval(tfit_rabi,tt);


%%% Uncertainty model
%   Constant error fraction
P_errfrac=harmmean(P_rabi_std(:,idx_mJ)./P_rabi_avg(:,idx_mJ));      % gives a reasonal fit to ~pi/2


%% DATA VISUALIZATION
%%% config
font_siz_reg=12;
font_siz_sml=10;
font_siz_lrg=14;
mark_siz=7;
line_wid=1.1;

[cc,clight,cdark]=palette(n_mf);
mark_typ={'o','^','d'};


%%% plot
figure('Name','rabi_oscillation');
hold on;

% Fitted model (a composite graphics object) - plotted first to be in bottom
lineProps.col={clight(idx_mJ,:)};
tfitted=mseb(1e6*tt,pp,P_errfrac.*pp,lineProps,1);

% DATA
h=NaN(n_mf,1);
for ii=1:n_mf
    th=ploterr(1e6*par_T,P_rabi_avg(:,ii),[],P_rabi_std(:,ii),'o','hhxy',0);
%     th=ploterr(tau_rot,P_rabi_avg(:,ii),[],P_rabi_std(:,ii),'o','hhxy',0);
    set(th(1),'color',cc(ii,:),'Marker',mark_typ{ii},'LineWidth',line_wid,...
        'MarkerSize',mark_siz,'MarkerFaceColor',clight(ii,:),...
        'DisplayName',num2str(configs.mf(ii).mf));
    set(th(2),'color',cc(ii,:),'LineWidth',line_wid);
    
    h(ii)=th(1);
end

% annotations
ax=gca;
set(ax,'FontSize',font_siz_reg);
set(ax,'Layer','Top');     % graphics axes should be always on top
axis tight;
box on;

ylim([0,1]);
ax.YTick=0:0.2:1;

xlabel('Pulse duration [$\mu$s]');
% xlabel('Pulse duration, $\tau$');
ylabel('$P$');


lgd=legend(h,'Location','East');
title(lgd,'$m_F$');
set(lgd,'FontSize',font_siz_reg);


%% Analysis for momentum resolved mJ-oscillation
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

%%% Cart-vecs to sph-polar vecs
%   * ks vectors are [th,phi,norm]
ks_par=cell(1,nparam);
for ii=1:nparam
    ks_par{ii} = cellfun(@(q) zxy2sphpol(q),k_par{ii},'UniformOutput',false);
end


%%% define momentum zones
% (A) lat-lon grid
nzone_th=4;     % num. zones to equipartition azimuthal (-pi,pi]
nzone_phi=2;    % for elevation [-pi/2,pi/2]
momzone_th=linspace(-pi,pi,nzone_th+1);
% momzone_phi=linspace(-pi/2,pi/2,nzone_phi+1);
momzone_phi=linspace(-0.8,0.8,nzone_phi+1);     % elev limits to BEC caps

% (B) try conical zone (either overlaps or misses)


%%% histogram atoms into zones
% (A) lat-lon grid: 2D histogram
hist_ed={momzone_th,momzone_phi};
nsc_zone=cell(1,nparam);
for ii=1:nparam
    nsc_zone{ii}=cell(1,n_mf);
    
    temp_nsc_zone=cellfun(@(ks) nhist(ks(:,1:2),hist_ed),ks_par{ii},'UniformOutput',false);
    
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


%%% DATA VIS
[ccc,ccclight,cccdark]=palette(nzone_th*nzone_phi);

h_rabi_momzone=figure('Name','rabi_momzone');
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
        set(th(1),'color',ccc(jj,:),'Marker',mark_typ{ii},'LineWidth',line_wid,...
            'MarkerSize',mark_siz,'MarkerFaceColor',ccclight(jj,:));      % 'DisplayName',num2str(configs.mf(ii).mf)
        set(th(2),'color',ccc(jj,:),'LineWidth',line_wid);
        
        h=cat(1,h,th(1));
    end
end

% annotations
ax=gca;
set(ax,'FontSize',font_siz_reg);
set(ax,'Layer','Top');     % graphics axes should be always on top
box on;

axis tight;
ax.YTick=0:0.2:1;
ylim([0,1]);

xlabel('Pulse duration [$\mu$s]');
% xlabel('Pulse duration, $\tau$');
ylabel('$P$');



%% Rabi oscillation (momentum-zone resolved fit)
momzone_amp=NaN(nzone_th,nzone_phi);
momzone_om=NaN(nzone_th,nzone_phi);

P_momzone_0=P_mJ_zone_avg{idx_mJ};  % get pop fracs for all (expparams,zones) for mJ=0

htemp=100;
for ii=1:numel(momzone_amp)
    [mm,nn]=ind2sub(size(momzone_amp),ii);  % get this zone
    tp=squeeze(P_momzone_0(mm,nn,:));    % pop fraction profile for this zone
    
    % check if this zone has any atoms
    if sum(isnan(tp))>0.33*length(tp)     % all NaN (a little too strict)
        continue
    end
        
    % estimate params
    tRabiAmp=max(tp)-min(tp);
    tRabiOmega=2*pi/(20e-5);        % we keep this near const
    
    tparam0=[tRabiAmp,tRabiOmega];
    
    
    %%% fit model
    tfopts=statset('Display','off');
    
    tfit_rabi=fitnlm(par_T,tp,rabi_mdl,tparam0,'CoefficientNames',{'amp','om'},...
        'Options',tfopts);
    
    tparam_fit=tfit_rabi.Coefficients.Estimate;
    
    % get fitted model params
    momzone_amp(mm,nn)=tparam_fit(1);
    momzone_om(mm,nn)=tparam_fit(2);
    
    % evaluate fitted model
    tt=linspace(0,max(par_T),1e3);
    pp=feval(tfit_rabi,tt);
    
    % vis
    figure(h_rabi_momzone);
    hold on;
%     plot(1e6*par_T,tp,'o','Color',ccc(ii,:));
    tfitted=plot(1e6*tt,pp,'-','Color',ccc(ii,:));
    uistack(tfitted,'bottom');    
end


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
% momzone_th=linspace(-pi,pi,nzone_th+1);
% momzone_phi=linspace(-pi/2,pi/2,nzone_phi+1);
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
    
    scatter_zxy(k_momzone_exp{mm,nn},plot_scat_size,ccc(ii,:));
    % color to match the zonal Rabi plot
end

view(3);
axis equal;
axis off;
