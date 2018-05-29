% Analysis of global rotation on correlated atom-pair 
%
% 2018.05.04: raw basis-dependent spin correlation - low mode occupancy
%
% DKS
%

%%% configs
path_src='C:\Users\HE BEC\bell\2018april\exp5_low_n_theta_v3\src';


%% load config
run(fullfile(path_src,'config_v1'));


%% load param log
if configs.flag.param_scan
    param_log=load_logfile(configs.path.paramlog);
    param_array = paramlog2array(param_log);
    
    % get unique param-vecs and tag each shot with param-ID
    [params,~,Ipar] = unique(param_array(:,2:end),'rows');
    
    param_id=param_array(:,1);
    nparam=size(params,1);      % number of unique param-set
    
    % group shot-ids by exp-param
    id_in_param=cell(1,nparam);
    for ii=1:nparam
        id_in_param{ii}=param_id(Ipar==ii);
    end
    
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

%% g2 BB
% MONEY

% TODO 
%   [ ] halo centering
%   [ ] improve filtering
%


g2=cell(nparam,1);
dk=cell(nparam,1);
for ii=1:nparam
    % run full g2
    [tg2,tdk]=summary_disthalo_g2(k_par{ii},0,true,0);     
    
    % store
    g2{ii}=tg2;
    dk{ii}=tdk;
end
    

%% Correlation coefficient
%    Note: may contain systematic error since unverified by Jan
%       + B-field distortions in halos aren't perfectly filtered
%


g2anti_par=NaN(nparam,1);
g2corr_par=NaN(nparam,1);
E_par=NaN(nparam,1);
E0_par=NaN(nparam,1);

for ii=1:nparam
    % get this paramset
    tg2=g2{ii};
    
    % aproximate g2 amplitude
%     g2corr_par(ii)=max([tg2{1}(15,15,15),tg2{2}(15,15,15)]);     
    g2corr_par(ii)=mean([tg2{1}(15,15,15),tg2{2}(15,15,15)]);     
    g2anti_par(ii)=tg2{3}(15,15,15);
    
    
    % evaluate spin-corr based on g2
    [E_par(ii),E0_par(ii)]=g2toE(g2corr_par(ii),g2anti_par(ii));
end



%% PRELIM bootstrapping
%%% CONFIG
subset_shotsize=300;    % shot-size of bootstrap sampled subset
n_subset=20;            % no. of bootstrap repeats
%   NOTE: unclear at the moment how config affects analysis


%%% main
nshot_par=shotSize(k_par);    % num. exp shots for each scanned parameter set
n_frac_samp=subset_shotsize./nshot_par; 


g2anti_samp=cell(nparam,1);
g2corr_samp=cell(nparam,1);
E_samp=cell(nparam,1);
E0_samp=cell(nparam,1);

for ii=1:nparam
    tnshot=nshot_par(ii);
    tn_frac_samp=n_frac_samp(ii);
    tk_par=k_par{ii};
    
    g2anti_samp{ii}=NaN(n_subset,1);
    g2corr_samp{ii}=NaN(n_subset,1);
    E_samp{ii}=NaN(n_subset,1);
    E0_samp{ii}=NaN(n_subset,1);
    
    for jj=1:n_subset
        Isamp=rand(tnshot,1)<tn_frac_samp;       % randomly select subset of shots to sample
        k_samp=tk_par(Isamp,:);           % the sampled data
        
        [tg2,tdk]=summary_disthalo_g2(k_samp,0,0,0);      % evaluate function
        
%         tg2corr=max([tg2{1}(15,15,15),tg2{2}(15,15,15)]);     % get results
        tg2corr=mean([tg2{1}(15,15,15),tg2{2}(15,15,15)]);     % get results
        tg2anti=tg2{3}(15,15,15);
        
        [tE,tE0]=g2toE(tg2corr,tg2anti);
        
        % store
        g2anti_samp{ii}(jj)=tg2anti;
        g2corr_samp{ii}(jj)=tg2corr;
        E_samp{ii}(jj)=tE;
        E0_samp{ii}(jj)=tE0;
    end
end

% statistics
E_bootstrap_mean=cellfun(@(x) mean(x,'omitnan'),E_samp);
E_bootstrap_sdev=cellfun(@(x) std(x,'omitnan'),E_samp);

E0_bootstrap_mean=cellfun(@(x) mean(x,'omitnan'),E0_samp);
E0_bootstrap_sdev=cellfun(@(x) std(x,'omitnan'),E0_samp);

%%% diplay some output
par_T'
E_par'
[E_bootstrap_mean,E_bootstrap_sdev]'
% [E0_bootstrap_mean,E0_bootstrap_sdev]'


%% Correlation - data vis
% Measured correlation
figure('Name','E_raw');

hh=ploterr(par_T/10e-6,E_par,...
    [],E_bootstrap_sdev,...
    'o','hhxy',0);
set(hh(1),'Marker','o','MarkerSize',7,...
    'Color','b','LineWidth',1.2,...
    'MarkerFaceColor','w');
set(hh(2),'Color','b','LineWidth',1.2);  % Y-err

% xlim([0,1]);
ylim([-1,1]);
xlabel('$\theta_0/\pi$');
ylabel('$E$');


% Mode occupancy corrected correlation
figure('Name','E_corrected');

hh=ploterr(par_T/10e-6,E0_par,...
    [],E0_bootstrap_sdev,...
    'o','hhxy',0);
set(hh(1),'Marker','o','MarkerSize',7,...
    'Color','k','LineWidth',1.2,...
    'MarkerFaceColor','w');
set(hh(2),'Color','k','LineWidth',1.2);  % Y-err

% xlim([0,1]);
ylim([-1,1]);
xlabel('$\theta_0/\pi$');
ylabel('$\bar{E}$');

ylim([-1,1]);


%% simple vis of 3D g2
%   NOTE: 3d g2 is completely averaged across Y-axis
for ii=1:nparam
    
    tT=1e6*par_T(ii);
    tfigname=sprintf('g2_T%0.3g',tT);
    
    figure('Name',tfigname);
    hold on;
    for jj=1:3      % 3-types of spin-spin type of mom-correlation
        subplot(1,3,jj);
        imagesc(mean(g2{ii}{jj},3));        % averaged across Y (3rd dim)
        
        title(sprintf('TH=%d, G2=%d',ii,jj));
        colorbar;
        
        axis square;
    end
end



%%
% %%% g2 spatial distribution
% % data to analyse
% k=k_par{1};
% 
% % config
% naz=200;
% nel=100;
% [az,el]=sphgrid(naz,nel);
% 
% b_pole=(abs(el)>asin(0.6));      % bool to bad region around poles
% 
% 
% % count
% NN=cellfun(@(x) haloZoneCount(x,az,el,0.03,[],'simple'),k,'UniformOutput',false);
% Nhalo{1}=cat(3,NN{:,1});
% Nhalo{2}=cat(3,NN{:,2});
% 
% Nhalo{1}(repmat(b_pole,[1,1,nshot]))=NaN;
% Nhalo{2}(repmat(b_pole,[1,1,nshot]))=NaN;
% 
% % analysis
% summary_disthalo_g2dist(Nhalo,az,el,1);
% 
% 
