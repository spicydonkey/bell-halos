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
% He*
C_gymag=2.8e6;     % gyromagnetic ratio (gamma) [Hz/G]

% vis
f_units='normalized';
f_pos=[0.2,0.2,0.2,0.3];
f_pos_wide=[0.2,0.2,0.25,0.3];
f_ren='painters';

cvir=viridis(3);
clvir=colshades(cvir);
[col,coll,cold]=palette(3);
c_gray=0.8*ones(1,3);
line_sty={'-','--',':','-.'};
mark_typ={'o','s','^','d'};
str_mm={'$x$','$y$','$z$'};
str_ss={'$\vert\!\uparrow\rangle$','$\vert\!\downarrow\rangle$'};
% str_ss={'$m_J = 1$','$m_J = 0$'};
% str_ss={'$\uparrow\uparrow$','$\downarrow\downarrow$','$\uparrow\downarrow$'};
mark_siz=7;
line_wid=1.5;
fontsize=12;
ax_lwidth=1.2;


%% load raw data
run('config_v1');


%% Get experimental params
% load logfile
param_log=load_logfile(configs.path.paramlog);
param_array = paramlog2array(param_log);
%   format: ID, T_PULSE, PHI, T_DELAY

% get unique param-vecs and tag each shot with param-ID
[params,~,Ipar] = unique(flip(param_array(:,2:end),2),'rows');
params=flip(params,2);
% sorted such that the Ith param (vec) can be binned into
%   T_PULSE X PHI X T_DELAY (3D) shaped array

param_id=param_array(:,1);
nparam=size(params,1);      % number of unique param-set

par_comp=arrayfun(@(c) unique(params(:,c)), 1:size(params,2),'UniformOutput',false);    % unique param-vec components
ncomp=cellfun(@(x) numel(x), par_comp);  % unique num vals for each comps in the set of param vecs
% ncomp=arrayfun(@(c) numel(unique(params(:,c))), 1:size(params,2));  % unique num vals for each comps in the set of param vecs

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


% DEBUG
% h_txy_raw=figure;
% plot_zxy(txy,1e5);
% axis equal;
% xlabel('x');
% ylabel('y');
% zlabel('t');
% view([0,0]);

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

% DEBUG
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
% k_par=cell(1,nparam);
k_par=cell(ncomp);
for ii=1:nparam
%     k_par{ii}=k_halo_filt(b_paramset(:,ii),:);      % get all halo data
    k_par{ii}=k_halo_filt(b_paramset{ii},:);      % get all halo data
    %from this param-set and store
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
%% Ramsey fringes
T=par_comp{3};      % T_delay between two pulses [s]
n_shot=cellfun(@(x) size(x,1),k_par,'UniformOutput',false);     % num shots each param

Nmf=cellfun(@(x) shotSize(x),k_par,'UniformOutput',false);      % num in mJ halos

% statistics
Nmf_avg=cellfun(@(n) mean(n,1),Nmf,'UniformOutput',false);      % avg num in mJ halo
Nmf_std=cellfun(@(n) std(n,0,1),Nmf,'UniformOutput',false);     % std num in mJ halos
Nmf_se=cellfun(@(s,N) s/sqrt(N),Nmf_std,n_shot,'UniformOutput',false);  % serr num in mJ halos

% pop fraction
p=cellfun(@(n) n/sum(n),Nmf_avg,'UniformOutput',false);         % pop fraction mJ=1/0


%% Model fit: Ramsey
P_ramsey=arrayfun(@(I) cat(1,p{1,I,:}),1:ncomp(2),'UniformOutput',false);

% config model
ramsey_mdl='y~c+amp*cos(om*x+phi)';
ramsey_cname={'amp','om','phi','c'};
ramsey_par0=[0.5,2*pi*1.5,0,0.5];

% fit
ramsey_fit=arrayfun(@(I) cellfun(@(P) fitnlm(1e6*T,P(:,I),ramsey_mdl,ramsey_par0,'CoefficientNames',ramsey_cname),P_ramsey,'UniformOutput',false),...
    1:2,'UniformOutput',false);
ramsey_fit=cat(1,ramsey_fit{:});
% format: MJ X PHI

% get params
ramsey_fpar=arrayfun(@(I) cellfun(@(f) f.Coefficients.Estimate(I),ramsey_fit),...
    1:numel(ramsey_cname),'UniformOutput',false);
ramsey_fpar=cat(3,ramsey_fpar{:});

ramsey_fparerr=arrayfun(@(I) cellfun(@(f) f.Coefficients.SE(I),ramsey_fit),...
    1:numel(ramsey_cname),'UniformOutput',false);
ramsey_fparerr=cat(3,ramsey_fparerr{:});
% MJ X PHI X PAR#

% get freq
om_ramsey=ramsey_fpar(:,:,2);
omerr_ramsey=ramsey_fparerr(:,:,2);

%% data vis
%% Phase & T_delay 
figname='two_pulse';
h=figure('Name',figname,'Units',f_units,'Position',[0.2,0.2,0.5,0.3],'Renderer',f_ren);
hold on;

for ii=1:ncomp(1)
    for jj=1:ncomp(2)
        subplot(ncomp(1),ncomp(2),sub2ind(ncomp([2,1]),jj,ii));
        tp=cat(1,p{ii,jj,:});
        
        hold on;
        for kk=1:2
            plot(1e6*T,tp(:,kk),...
                'Color',clvir(kk,:),'LineStyle',line_sty{kk},...   %'LineWidth',line_wid,...
                'Marker',mark_typ{kk},'MarkerEdgeColor',cvir(kk,:),'MarkerFaceColor',clvir(kk,:),...       %'MarkerSize',mark_siz,...
                'DisplayName',str_ss{kk});
        end
        
        % annotate subplot
        ax=gca;
        box on;
        ax_ylim=ax.YLim;
        ax_xlim=ax.XLim;
        
        set(ax,'Layer','Top');
        ax.FontSize=fontsize;
        ax.LineWidth=ax_lwidth;
        xlabel('$T~[\mu s]$');
        ylabel('$P$');
    end
end

%% Ramsey fringe
figname='ramsey_fringe';
h=figure('Name',figname,'Units',f_units,'Position',[0.2,0.2,0.5,0.2],'Renderer',f_ren);
hold on;

xx=1e6*linspace(min(T),max(T));     % x-axis range for fitted curve

for ii=1:ncomp(2)
    subplot(1,ncomp(2),ii);
    tp=P_ramsey{ii};
    
    hold on;
    for jj=1:2
        pexp=plot(1e6*T,tp(:,jj),...
            'Color',clvir(jj,:),'LineStyle','none',...   %'LineWidth',line_wid,...
            'Marker',mark_typ{jj},'MarkerEdgeColor',cvir(jj,:),'MarkerFaceColor',clvir(jj,:),...       %'MarkerSize',mark_siz,...
            'DisplayName',str_ss{jj});
        
        yy=feval(ramsey_fit{jj,ii},xx);
        pfit=plot(xx,yy,'LineStyle',line_sty{jj},'Color',clvir(jj,:));
        uistack(pfit,'bottom');
    end
    titlestr=sprintf('%s %0.3g(%0.1g) MHz','$f_L=$',om_ramsey(1,ii)/(2*pi),omerr_ramsey(1,ii)/(2*pi));
    title(titlestr);
    
    % annotate subplot
    ax=gca;
    box on;
    ax_ylim=ax.YLim;
    ax_xlim=ax.XLim;
    
    set(ax,'Layer','Top');
    ax.FontSize=fontsize;
    ax.LineWidth=ax_lwidth;
    xlabel('Pulse delay $T~[\mu s]$');
    ylabel('Pop fraction $P$');
    xlim(1e6*[min(T),max(T)]+[-0.1,0.1]);
    ylim([-0.05,1.05]);
end