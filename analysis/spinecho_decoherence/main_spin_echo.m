% spin-echo experiment
%   spin-echo for ~long delay
%   
%   approx optimal polarization to eliminate sigma+ component
%   global rotation - beam very large compared to spatial distribution of atoms
%   rotation on (1,1) halo source
%
%   
%
%   NOTE:
%       * significant instability in PUSH/NULLER-SO sequence: mf=1 position variability
%
% 20180409
% DKS
%


%% load raw data
run('config_v1');


%% Get experimental params
% load wfmgen log
[params,id_in_param,param_id,Ipar]=wfmgen_log_parser(configs.path.paramlog);
nparam=size(params,1);      % number of unique param-set

% get searched param
Tdelay=params;


%% Load txy
[txy,fout]=load_txy(configs.load.path,configs.load.id,configs.load.window,configs.load.mincount,configs.load.maxcount,[],1,0,0);

zxy=txy2zxy(txy,vz);
shot_id=fout.id_ok;

% build boolean array selector for scanned param-set
b_paramset=cellfun(@(idx) ismember(shot_id,idx),...
    id_in_param,'UniformOutput',false);
b_paramset=horzcat(b_paramset{:});


% DEBUG
h_txy_raw=figure;
plot_zxy(txy,1e5);
axis equal;
xlabel('x');
ylabel('y');
zlabel('t');
view([0,0]);

% h_zxy_raw=figure;
% plot_zxy(zxy,1e5);
% axis equal;
% xlabel('x');
% ylabel('y');
% zlabel('z');
% view([0,0]);


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
h_zxy_raw=figure();
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
clear txy zxy zxy0 zxy0_filt tzxy tzxy_mf tp_bec tp_bec0 tp_halo tr_bec0 tr_lim;
clear h_zxy*;       % clear figs


%% ANALYSIS
% non-dimensionalised time: Larmor period
T_larmor=0.65e-6;       % larmor precession period (s)
dtau=Tdelay/T_larmor;


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


% %% Number counting
% Nmf=cellfun(@(x) shotSize(x),k_par,'UniformOutput',false);
% 
% % statistics
% Nmf_avg=cellfun(@(n) mean(n,1),Nmf,'UniformOutput',false);
% Nmf_avg=vertcat(Nmf_avg{:});
% 
% Nmf_std=cellfun(@(n) std(n,0,1),Nmf,'UniformOutput',false);
% Nmf_std=vertcat(Nmf_std{:});
% 
% Nmf_se=Nmf_std./sqrt(cellfun(@(x)size(x,1),Nmf))';


% % population fraction
% Nmf_tot=sum(Nmf_avg,2);         % total number in all spin-state
% 
% P=Nmf_avg./Nmf_tot;
% P_err=(vnorm(Nmf_se,2)./Nmf_tot);      % simple estimate of error
% 
% % vis
% figure('Name','population_fraction');
% hold on;
% 
% p1=ploterr(dtau,P(:,1),[],P_err,'ro');
% p2=ploterr(dtau,P(:,2),[],P_err,'bo');
% 
% p1(1).DisplayName='1';
% p2(1).DisplayName='0';
% 
% titlestr=sprintf('spin-echo P-fraction');
% title(titlestr);
% xlabel('$\tau$');
% ylabel('P');
% box on;
% ax=gca;
% xlim([0,ax.XLim(2)]);
% ylim([0,ax.YLim(2)]);
% 
% lgd=legend([p1(1),p2(1)]);
% lgd.Title.String='$m_F$';


%% DATA VISUALIZATION
% config
font_siz_reg=12;
font_siz_sml=10;
font_siz_lrg=14;
mark_siz=7;
line_wid=1.1;
[c0,clight,cdark]=palette(n_mf);
mark_typ={'o','^','d'};

%%% plot
hf=figure('Name','spinecho_halo');
hold on;


% data
h=NaN(n_mf,1);
for ii=1:n_mf
    th=ploterr(1e6*Tdelay,P_mJ_halo_avg(:,ii),[],P_mJ_halo_std(:,ii),'o','hhxy',0);
%     th=ploterr(tau_rot,P_mJ_halo_avg(:,ii),[],P_mJ_halo_std(:,ii),'o','hhxy',0);
    set(th(1),'color',c0(ii,:),'Marker',mark_typ{ii},'LineWidth',line_wid,...
        'MarkerSize',mark_siz,'MarkerFaceColor',clight(ii,:),...
        'DisplayName',num2str(configs.mf(ii).mf));
    set(th(2),'color',c0(ii,:),'LineWidth',line_wid);
    
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

xlabel('Pulse delay [$\mu$s]');
% xlabel('Pulse delay, $\tau$');
ylabel('$P$');

lgd=legend(h,'Location','East');
title(lgd,'$m_F$');
set(lgd,'FontSize',font_siz_reg);


%% Model: spin-echo and decay
%   TODO
% pp=P_mJ_halo_avg(:,1);
% 
% % v1
% modelfun = @(b,x) b(1)*exp(-b(2)*x);
% beta0 = [1,1e-3];
% 
% % v2
% % modelfun = @(b,x) 1*exp(-b(1)*x);
% % beta0 = [1e-3];
% 
% mdl = fitnlm(dtau,pp,modelfun,beta0);
% 
% %%% eval fit
% beta_fit=mdl.Coefficients.Estimate;
% 
% tfit=linspace(0,700,1e5);
% pfit=modelfun(beta_fit,tfit);
% 
% % % estimate
% % beta_est=[0.967,0.00015];
% % pest=modelfun(beta_est,tfit);
% 
% %%% vis
% figure('Name','population_fraction');
% hold on;
% 
% p1=ploterr(dtau,P_mJ_halo_avg(:,1),[],P_mJ_halo_std(:,1),'ro');
% p1(1).DisplayName='data';
% 
% titlestr=sprintf('Spin-echo');
% title(titlestr);
% xlabel('$\tau$');
% ylabel('$P_1$');
% box on;
% % ax=gca;
% % xlim([0,ax.XLim(2)]);
% % ylim([0,ax.YLim(2)]);
% 
% p_mdl=plot(tfit,pfit,'k-','DisplayName','fit');
% uistack(p_mdl,'bottom')
% 
% % p_est=plot(tfit,pest,'--','Color',0.3*[1,1,1],'DisplayName','estimate');
% % uistack(p_mdl,'bottom')
% 
% lgd=legend([p1(1),p_mdl]);
% % lgd.Title.String='$m_F$';
% 
% ax=gca;
% xlim([0,ax.XLim(2)]);
% ylim([0,ax.YLim(2)]);


%% MOMENTUM-RESOLVED: mJ-oscillation
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


%% DATA VIS
[cc0,cclight,ccdark]=palette(nzone_th*nzone_phi);

h_spinecho_momzone=figure('Name','spinecho_momzone');
hold on;

h=[];
for ii=1:n_mf
    % iterate over lat-long grid 2D indices
    for jj=1:nzone_th*nzone_phi
        [mm,nn]=ind2sub([nzone_th,nzone_phi],jj);
        th=ploterr(1e6*Tdelay,squeeze(P_mJ_zone_avg{ii}(mm,nn,:)),[],...
            squeeze(P_mJ_zone_std{ii}(mm,nn,:)),'o','hhxy',0);
%         th=ploterr(dtau,squeeze(P_mJ_zone_avg{ii}(mm,nn,:)),[],...
%             squeeze(P_mJ_zone_std{ii}(mm,nn,:)),'o','hhxy',0);
        set(th(1),'color',cc0(jj,:),'Marker',mark_typ{ii},'LineWidth',line_wid,...
            'MarkerSize',mark_siz,'MarkerFaceColor',cclight(jj,:));      % 'DisplayName',num2str(configs.mf(ii).mf)
        set(th(2),'color',cc0(jj,:),'LineWidth',line_wid);
        
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

xlabel('Pulse delay [$\mu$s]');
% xlabel('Pulse delay, $\tau$');
ylabel('$P$');



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
    
    scatter_zxy(k_momzone_exp{mm,nn},plot_scat_size,cc0(ii,:));
    % color to match the zonal Rabi plot
end

view(3);
axis equal;
axis off;
