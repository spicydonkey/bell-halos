% Ramsey interferometry with L3 Raman beam
%
%
% DKS
% 2018-06-06


config_name='C:\Users\HE BEC\Documents\MATLAB\bell-halos\analysis\expS4_ramsey\src\config_1.m';


%% load config
run(config_name);


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
    par_dphi=params;       % scanned relative phase
else
    % TODO
    %   do I need to set some things to default? or re-code analysis?
    
    % defaults autoscan-related vars
    nparam=1;       % 1- since we don't search any params
    id_in_param={configs.load.id};    % all IDS to load
    
    par_dphi=NaN;      % default param to NaN
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


%% Fit Rabi oscillation
% TODO


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
figure('Name','ramsey_fringe');
hold on;


h=NaN(n_mf,1);
for ii=1:n_mf
    th=ploterr(par_dphi,P_rabi_avg(:,ii),[],P_rabi_std(:,ii),'o','hhxy',0);
    set(th(1),'color',cc(ii,:),'Marker',mark_typ{ii},'LineWidth',line_wid,...
        'MarkerSize',mark_siz,'MarkerFaceColor',clight(ii,:),...
        'DisplayName',num2str(configs.mf(ii).mf));
    set(th(2),'color',cc(ii,:),'LineWidth',line_wid);
    
    h(ii)=th(1);
end

set(gca,'FontSize',font_siz_reg);

set(gca,'Layer','Top');     % graphics axes should be always on top
box on;

xlabel('$\phi$');
% ylabel('$P\left(m_F\right)$');
ylabel('$P$');


lgd=legend(h,'Location','East');
title(lgd,'$m_F$');
set(lgd,'FontSize',font_siz_reg);

xlim([0-2*pi/25,2*pi+2*pi/25]);

% X-ticks
xticks(0:pi/2:2*pi);
xticklabels({'$0$','$\pi/2$','$\pi$','$3\pi/2$','$2\pi$'});