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
% save
do_save_figs=false;
dir_save='C:\Users\HE BEC\Dropbox\PhD\projects\halo_metrology\analysis\ramsey\prelim_20181128';

% He*
C_gymag=2.8e6;     % gyromagnetic ratio (gamma) [Hz/G]

% vis
f_units='normalized';
f_pos=[0.2,0.2,0.2,0.3];
f_pos_wide=[0.2,0.2,0.25,0.3];
f_ren='painters';

cvir=viridis(3);
clvir=colshades(cvir);
cvir2=viridis(4);
clvir2=colshades(cvir2);
[col,coll,cold]=palette(3);
c_gray=0.8*ones(1,3);
line_sty={'-','--',':','-.'};
mark_typ={'o','s','^','d'};
% str_mm={'$x$','$y$','$z$'};
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
phi=par_comp{2};    % phase delay of 2nd pi/2 pulse (rad)

n_shot=cellfun(@(x) size(x,1),k_par,'UniformOutput',false);     % num shots each param

%%% ATOM NUMBERS
Nm_shot=cellfun(@(x) shotSize(x),k_par,'UniformOutput',false);      % num in mJ halos

% statistics
Nm_avg=cellfun(@(n) mean(n,1),Nm_shot,'UniformOutput',false);      % avg num in mJ halo
Nm_std=cellfun(@(n) std(n,0,1),Nm_shot,'UniformOutput',false);     % std num in mJ halos
Nm_se=cellfun(@(s,N) s/sqrt(N),Nm_std,n_shot,'UniformOutput',false);  % serr num in mJ halos

N_tot=cellfun(@(n) sum(n),Nm_avg);     % avg total atom num det'd (mJ summed)
N_tot_avg=mean(N_tot,3);                % avg total atom num det'd (T avg'd)

%%% pop fraction (pop frac per spin comp)
pm_shot=cellfun(@(n) n./sum(n,2),Nm_shot,'UniformOutput',false);
pm_avg=cellfun(@(x) mean(x,1),pm_shot,'UniformOutput',false);
pm_std=cellfun(@(x) std(x,0,1),pm_shot,'UniformOutput',false);
pm_se=cellfun(@(s,N) s/sqrt(N),pm_std,n_shot,'UniformOutput',false);

%%% Population (inversion): P = p_up - p_down
P_shot=cellfun(@(n) (n(:,1)-n(:,2))./sum(n,2),Nm_shot,'UniformOutput',false);   
P_avg=cellfun(@(x) mean(x,1),P_shot,'UniformOutput',false);
P_std=cellfun(@(x) std(x,0,1),P_shot,'UniformOutput',false);
P_se=cellfun(@(s,N) s/sqrt(N),P_std,n_shot,'UniformOutput',false);

% TODO 20190110
P_ramsey=arrayfun(@(I) cat(1,P_avg{1,I,:}),1:ncomp(2),'UniformOutput',false);
P_std_ramsey=arrayfun(@(I) cat(1,P_std{1,I,:}),1:ncomp(2),'UniformOutput',false);
P_err_ramsey=arrayfun(@(I) cat(1,P_se{1,I,:}),1:ncomp(2),'UniformOutput',false);

figname='2_pi/2_pulse_with_delay';
h=figure('Name',figname,'Units','centimeters','Position',[0,0,8.6,4],'Renderer','Painters');
hold on;
for ii=1:ncomp(2)
    plot(1e6*T,P_ramsey{ii},'o-');
end
box on;
xlabel('$\tau$ ($\mu$s)');
ylabel('P');
lgd=legend(num2str(phi/pi),'Location','East');
title(lgd,'$\phi/\pi$');

%% Model fit: Ramsey
%%% pop fraction: p±
p_ramsey=arrayfun(@(I) cat(1,pm_avg{1,I,:}),1:ncomp(2),'UniformOutput',false);
p_std_ramsey=arrayfun(@(I) cat(1,pm_std{1,I,:}),1:ncomp(2),'UniformOutput',false);
p_err_ramsey=arrayfun(@(I) cat(1,pm_se{1,I,:}),1:ncomp(2),'UniformOutput',false);

% config model
ramsey_mdl='y~c+amp*cos(om*x+phi)';
ramsey_cname={'amp','om','phi','c'};
ramsey_par0=[0.5,2*pi*1.5,0,0.5];

% fit
ramsey_fit=arrayfun(@(I) cellfun(@(P) fitnlm(1e6*T,P(:,I),ramsey_mdl,ramsey_par0,'CoefficientNames',ramsey_cname),p_ramsey,'UniformOutput',false),...
    1:2,'UniformOutput',false);
ramsey_fit=cat(1,ramsey_fit{:});
% format: MJ X PHI_DELAY

% get params
ramsey_fpar=arrayfun(@(I) cellfun(@(f) f.Coefficients.Estimate(I),ramsey_fit),...
    1:numel(ramsey_cname),'UniformOutput',false);
ramsey_fpar=cat(3,ramsey_fpar{:});

ramsey_fparerr=arrayfun(@(I) cellfun(@(f) f.Coefficients.SE(I),ramsey_fit),...
    1:numel(ramsey_cname),'UniformOutput',false);
ramsey_fparerr=cat(3,ramsey_fparerr{:});
% MJ X PHI_DELAY X PAR#

% get freq
om_ramsey=ramsey_fpar(:,:,2);
omerr_ramsey=ramsey_fparerr(:,:,2);

% magnetic field
B=1e6*om_ramsey/(2*pi*C_gymag);
Berr=1e6*omerr_ramsey/(2*pi*C_gymag);


%% VIS: Phase & T_delay 
figname='two_pulse';
h=figure('Name',figname,'Units',f_units,'Position',[0.2,0.2,0.5,0.3],'Renderer',f_ren);
hold on;

for ii=1:ncomp(1)
    for jj=1:ncomp(2)
        subplot(ncomp(1),ncomp(2),sub2ind(ncomp([2,1]),jj,ii));
        tp=cat(1,pm_avg{ii,jj,:});
        
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
        
        xlim(1e6*[min(T),max(T)]+[-0.1,0.1]);
        ylim([-0.05,1.05]);
    end
end

% save fig
if do_save_figs
    savefigname=sprintf('fig_%s_%s',figname,getdatetimestr);
    fpath=fullfile(dir_save,savefigname);
    
    saveas(h,strcat(fpath,'.fig'),'fig');
    print(h,strcat(fpath,'.svg'),'-dsvg');
end

%% VIS: Ramsey fringe
figname='ramsey_fringe';
h=figure('Name',figname,'Units',f_units,'Position',[0.2,0.2,0.5,0.2],'Renderer',f_ren);
hold on;

xx=1e6*linspace(min(T),max(T),1e3);     % x-axis range for fitted curve

for ii=1:ncomp(2)
    subplot(1,ncomp(2),ii);
    tp=p_ramsey{ii};
%     tperr=Pstd_ramsey{ii};
    tperr=p_err_ramsey{ii};
    
    hold on;
    for jj=1:2        
        pexp=ploterr(1e6*T,tp(:,jj),[],tperr(:,jj),mark_typ{jj},'hhxy',0);
        set(pexp(1),'MarkerFaceColor',clvir(jj,:),'MarkerEdgeColor',cvir(jj,:),...
            'DisplayName',str_ss{jj});
        set(pexp(2),'Color',cvir(jj,:));
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

% save fig
if do_save_figs
    savefigname=sprintf('fig_%s_%s',figname,getdatetimestr);
    fpath=fullfile(dir_save,savefigname);
    
    saveas(h,strcat(fpath,'.fig'),'fig');
    print(h,strcat(fpath,'.svg'),'-dsvg');
end

%% MODE RESOLVED RAMSEY
%%% construct spatial zones at latlon grid + solid angle
alpha=pi/10;         % half-cone angle
lim_az=[-pi,pi];    % no inversion symmetry
phi_max=pi/4;       
lim_el=[-phi_max,phi_max];

n_az=20;                	% equispaced bins
n_el=5;

az_disp=deg2rad(-180:90:90);     % azim sections (great circles) to display
% el_disp=deg2rad(-30:30:30);    % elev/lat zones to display

az=linspace(lim_az(1),lim_az(2),n_az+1);
az=az(1:end-1);
el=linspace(lim_el(1),lim_el(2),n_el);

[~,iaz_disp]=arrayfun(@(x) min(abs(az-x)),az_disp);     % idx to displayable azim
naz_disp=length(iaz_disp);

[gaz,gel]=ndgrid(az,el);    % az-el grid
n_zone=numel(gaz);

[~,iel_0]=min(abs(el));        % idx to ~zero elev angle (equator)

%%% get atom numbers in regions
k_ramsey=squeeze(k_par(1,:,:));     % halo k-vectors Ramsey exp
N_mode=cell(ncomp(2),ncomp(3));     % num atoms in zone categorise by exp param

for ii=1:numel(k_ramsey)
    tk=k_ramsey{ii};        % exp data for this param
    [iph,iT]=ind2sub([npar_phidelay,npar_tdelay],ii);
    
    % get num in zone/shot
    tN=arrayfun(@(th,ph) cellfun(@(x) size(inCone(x,th,ph,alpha),1),tk),...
        gaz,gel,'UniformOutput',false);     
    
    % format into multi-dim array: AZ X EL X SHOT X MJ
    ttN=NaN(cat(2,size(gaz),size(tN{1})));
    for jj=1:n_zone
        [iaz,iel]=ind2sub([n_az,n_el],jj);
        
        ttN(iaz,iel,:,:)=tN{iaz,iel};
    end
    N_mode{ii}=ttN;
end
% clearvars tN ttN;

% pop fraction
p_shot_mode=cellfun(@(x) x./sum(x,4),N_mode,'UniformOutput',false);
p_avg_mode=cellfun(@(x) mean(x,3),p_shot_mode,'UniformOutput',false);
p_std_mode=cellfun(@(x) std(x,0,3),p_shot_mode,'UniformOutput',false);
n_shot_mode=cellfun(@(x) size(x,3), p_shot_mode);

% form into multidim array: PHI_DELAY X T_DELAY X AZ X EL X [] X MJ
p_avg_mode=cell2mat(cellfun(@(a) reshape(a,[1,1,size(a)]),p_avg_mode,'UniformOutput',false));
p_std_mode=cell2mat(cellfun(@(a) reshape(a,[1,1,size(a)]),p_std_mode,'UniformOutput',false));

p_se_mode=p_std_mode./sqrt(n_shot_mode);        % standard error

%% analyse Ramsey fringe for all modes and exp params
ramsey_fit_mode=cell(npar_phidelay,n_az,n_el,2);
% format: PHI_DELAY X AZ X EL X MJ
for ii=1:npar_phidelay
    for jj=1:n_zone
        [iaz,iel]=ind2sub([n_az,n_el],jj);
        for kk=1:2
            tP=p_avg_mode(ii,:,iaz,iel,1,kk);       % pop oscillation
            ramsey_fit_mode{ii,iaz,iel,kk}=fitnlm(1e6*T,tP,ramsey_mdl,...
                ramsey_par0,'CoefficientNames',ramsey_cname);
        end
    end
end

% get fit params
ramsey_fpar_mode=arrayfun(@(I) cellfun(@(f) f.Coefficients.Estimate(I),ramsey_fit_mode),...
    1:numel(ramsey_cname),'UniformOutput',false);
ramsey_fpar_mode=cat(ndims(ramsey_fpar_mode{1})+1,ramsey_fpar_mode{:});

ramsey_fparerr_mode=arrayfun(@(I) cellfun(@(f) f.Coefficients.SE(I),ramsey_fit_mode),...
    1:numel(ramsey_cname),'UniformOutput',false);
ramsey_fparerr_mode=cat(ndims(ramsey_fparerr_mode{1})+1,ramsey_fparerr_mode{:});
% PHI_DELAY X AZ X EL X MJ X PAR#

% get freq
om_ramsey_mode=ramsey_fpar_mode(:,:,:,:,2);
omerr_ramsey_mode=ramsey_fparerr_mode(:,:,:,:,2);
% PHI_DELAY X AZ X EL X MJ

% outlier to NaN
b_outlier=isoutlier(om_ramsey_mode);    % find outlier by Med Abs Dev
om_ramsey_mode_raw=om_ramsey_mode;      % store original raw data set
omerr_ramsey_mode_raw=omerr_ramsey_mode;
om_ramsey_mode(b_outlier)=NaN;          % outlier --> NaN
omerr_ramsey_mode(b_outlier)=NaN;

om_mode_0=squeeze(mean(om_ramsey_mode,1,'omitnan'));   % mode-resolved L-freq avgd thru PHI_DELAY
n_samp=squeeze(sum(~b_outlier,1));       % num exp samples to average over
omerr_mode_0=sqrt(squeeze(sum(omerr_ramsey_mode.^2,1,'omitnan')))./n_samp;

% magnetic field
B_mode=1e6*om_ramsey_mode/(2*pi*C_gymag);
Berr_mode=1e6*omerr_ramsey_mode/(2*pi*C_gymag);

B_mode_0=1e6*om_mode_0/(2*pi*C_gymag);          % average PHI_DELAY exp params
Berr_mode_0=1e6*omerr_mode_0/(2*pi*C_gymag);


%% VIS: Mode-resolved Ramsey fringe
figname='ramsey_fringe_mode';
h=figure('Name',figname,'Units',f_units,'Position',[0.2,0.2,0.5,0.2],'Renderer',f_ren);
hold on;

xx=1e6*linspace(min(T),max(T),1e3);     % x-axis range for fitted curve

idx_phi_disp=1;

for ii=1:naz_disp
    iaz=iaz_disp(ii);
    tp=squeeze(p_avg_mode(idx_phi_disp,:,iaz,iel_0,1,:));   % T X MJ
    tperr=squeeze(p_se_mode(idx_phi_disp,:,iaz,iel_0,1,:));
    
    subplot(1,naz_disp,ii);
    hold on;
    for jj=1:2
        pexp=ploterr(1e6*T,tp(:,jj),[],tperr(:,jj),mark_typ{jj},'hhxy',0);
        set(pexp(1),'MarkerFaceColor',clvir(jj,:),'MarkerEdgeColor',cvir(jj,:),...
            'DisplayName',str_ss{jj});
        set(pexp(2),'Color',cvir(jj,:));
        
        yy=feval(ramsey_fit_mode{idx_phi_disp,iaz,iel_0,jj},xx);
        pfit=plot(xx,yy,'LineStyle',line_sty{jj},'Color',clvir(jj,:));
        uistack(pfit,'bottom');
    end
    titlestr=sprintf('%s %0.0f','$\theta=$',rad2deg(az(iaz)));
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

% save fig
if do_save_figs
    savefigname=sprintf('fig_%s_%s',figname,getdatetimestr);
    fpath=fullfile(dir_save,savefigname);
    
    saveas(h,strcat(fpath,'.fig'),'fig');
    print(h,strcat(fpath,'.svg'),'-dsvg');
end

%% VIS: Magnetometry
idx_mJ=1;       % just use fit param from mJ=1

figname='halo_magnetometry';
h=figure('Name',figname,'Units',f_units,'Position',[0.2,0.2,0.35,0.5],'Renderer',f_ren);

for ii=1:npar_phidelay
    subplot(2,2,ii);
    
    tp=plotFlatMapWrappedRad(gaz,gel,squeeze(B_mode(ii,:,:,idx_mJ)),'eckert4','texturemap');
    
    % annotation
    ax=gca;
    set(ax,'Layer','Top');
    % ax.FontSize=fontsize;
    % ax.LineWidth=ax_lwidth;
    
    cbar=colorbar('SouthOutside');
    cbar.TickLabelInterpreter='latex';
    cbar.Label.Interpreter='latex';
    cbar.Label.String='Magnetic field $\mathrm{B}$ [G]';
    cbar.Label.FontSize=fontsize;
    cbar.FontSize=fontsize;
    colormap('viridis');
end

% save fig
if do_save_figs
    savefigname=sprintf('fig_%s_%s',figname,getdatetimestr);
    fpath=fullfile(dir_save,savefigname);
    
    saveas(h,strcat(fpath,'.fig'),'fig');
    print(h,strcat(fpath,'.png'),'-dpng','-r300');
end


%% VIS: Magnetic tomography: equatorial
figname='B_equatorial_tomography';
h=figure('Name',figname,'Units',f_units,'Position',f_pos,'Renderer',f_ren);

B_eq=B_mode(:,:,iel_0,idx_mJ);      % magnetic field around equator
Berr_eq=Berr_mode(:,:,iel_0,idx_mJ);      

hold on;
pleg=NaN(npar_phidelay,1);
for ii=1:npar_phidelay
    tp=ploterr(rad2deg(az),B_eq(ii,:),[],Berr_eq(ii,:),'o','hhxy',0);
    set(tp(1),'Marker',mark_typ{ii},...
        'MarkerFaceColor',clvir2(ii,:),'MarkerEdgeColor',cvir2(ii,:),...
        'DisplayName',num2str(rad2deg(par_comp{2}(ii)),3));
    set(tp(2),'Color',cvir2(ii,:));
    pleg(ii)=tp(1);
end

% annotation
box on;
ax=gca;
set(ax,'Layer','Top');
ax.FontSize=fontsize;
ax.LineWidth=ax_lwidth;

title('Equatorial tomography $\phi=0$');

xlabel('Azimuthal angle $\theta$ [$^\circ$]');
ylabel('Magnetic field $\mathrm{B}$ [G]');

xlim([-180,180]);
ax.XTick=-180:90:180;
% ylim([0.52,0.545]);

% lgd=legend(pleg);

% save fig
if do_save_figs
    savefigname=sprintf('fig_%s_%s',figname,getdatetimestr);
    fpath=fullfile(dir_save,savefigname);
    
    saveas(h,strcat(fpath,'.fig'),'fig');
    print(h,strcat(fpath,'.svg'),'-dsvg');
end


%% VIS: Magnetometry (PHI-avg)
figname='halo_magnetometry_phi0';
h=figure('Name',figname,'Units',f_units,'Position',f_pos,'Renderer',f_ren);

tp=plotFlatMapWrappedRad(gaz,gel,B_mode_0(:,:,idx_mJ),'eckert4','texturemap');

% annotation
ax=gca;
set(ax,'Layer','Top');
% ax.FontSize=fontsize;
% ax.LineWidth=ax_lwidth;

cbar=colorbar('SouthOutside');
cbar.TickLabelInterpreter='latex';
cbar.Label.Interpreter='latex';
cbar.Label.String='Magnetic field $\mathrm{B}$ [G]';
cbar.Label.FontSize=fontsize;
cbar.FontSize=fontsize;
colormap('viridis');

% save fig
if do_save_figs
    savefigname=sprintf('fig_%s_%s',figname,getdatetimestr);
    fpath=fullfile(dir_save,savefigname);
    
    saveas(h,strcat(fpath,'.fig'),'fig');
    print(h,strcat(fpath,'.png'),'-dpng','-r300');
end


%% VIS: Magnetic tomography: equatorial (PHI-avg)
figname='B_equatorial_tomography_phi0';
h=figure('Name',figname,'Units',f_units,'Position',f_pos,'Renderer',f_ren);

B_eq_0=B_mode_0(:,iel_0,idx_mJ);      % magnetic field around equator
Berr_eq_0=Berr_mode_0(:,iel_0,idx_mJ);      

hold on;
tp=ploterr(rad2deg(az),B_eq_0,[],Berr_eq_0,'o','hhxy',0);
set(tp(1),'Marker',mark_typ{2},...
    'MarkerFaceColor',clvir2(2,:),'MarkerEdgeColor',cvir2(2,:),...
    'DisplayName',num2str(rad2deg(par_comp{2}(2)),3));
set(tp(2),'Color',cvir2(2,:));
pleg=tp(1);


% annotation
box on;
ax=gca;
set(ax,'Layer','Top');
ax.FontSize=fontsize;
ax.LineWidth=ax_lwidth;

title('Equatorial tomography $\phi=0$');

xlabel('Azimuthal angle $\theta$ [$^\circ$]');
ylabel('Magnetic field $\mathrm{B}$ [G]');

% axis tight;
xlim([-180,180]);
xlim([-180,180]);
ax.XTick=-180:90:180;
% ylim([0.52,0.545]);

% lgd=legend(pleg);

% save fig
if do_save_figs
    savefigname=sprintf('fig_%s_%s',figname,getdatetimestr);
    fpath=fullfile(dir_save,savefigname);
    
    saveas(h,strcat(fpath,'.fig'),'fig');
    print(h,strcat(fpath,'.svg'),'-dsvg');
end