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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau=par_comp{3};              % T_delay between two pulses [s]
phi=par_comp{2};            % phase delay of 2nd pi/2 pulse (rad)

% Pre-processing
idx_ramsey_exp=1;           % 1: Ramsey (pi/2), 2: pi-pulses
idx_phi0=find(phi==0);      % two-pi/2-pulse (Ramsey) dataset with zero phase delay phi=0

k_ramsey_phi0=squeeze(k_par(idx_ramsey_exp,idx_phi0,:));     % Ramsey sub-dataset (phi=0)
n_shot=cellfun(@(x) size(x,1),k_ramsey_phi0);       % number of shots per tau


%% atom number and population
% raw #spins: #atom in mJ halos
Nm_halo=cellfun(@(x) shotSize(x),k_ramsey_phi0,'UniformOutput',false);      

% #spin in halo
Nm_halo_avg=cellfun(@(n) mean(n,1),Nm_halo,'UniformOutput',false);      % avg num in mJ halo
Nm_halo_avg=vertcat(Nm_halo_avg{:});        % form as array
Nm_halo_std=cellfun(@(n) std(n,0,1),Nm_halo,'UniformOutput',false);     % std 
Nm_halo_std=vertcat(Nm_halo_std{:});
Nm_halo_se=Nm_halo_std./sqrt(n_shot);       % std-err

Ninv_halo=cellfun(@(x) x(:,1)-x(:,2),Nm_halo,'UniformOutput',false);    % "collective spin"/inversion

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


%% Ramsey model
%%% models
% pop fraction per spin comp
mramsey_mdl='y~c+amp*cos(om*x+phi)';
mramsey_cname={'amp','om','phi','c'};
mramsey_par0=[0.5,2*pi*1.5,0,0.5];

idx_om_mramsey=find(cellfun(@(s) strcmp(s,'om'),mramsey_cname)==1);     % param-vec idx to 'om'

% pop-inversion
Pramsey_mdl='y~amp*cos(om*x+phi)';
Pramsey_cname={'amp','om','phi'};
Pramsey_par0=[0.5,2*pi*1.5,0];

idx_om_Pramsey=find(cellfun(@(s) strcmp(s,'om'),Pramsey_cname)==1);


%% fit to halo-avgd data
%%% pop fraction per spin comp
% collate data to fit
tt=arrayfun(@(x,n) x*ones(n,1),tau,n_shot,'UniformOutput',false);
tt=vertcat(tt{:});
yy=vertcat(pm_halo{:});

% fit
mramsey_fit_halo=cell(1,2);
for ii=1:2
    mramsey_fit_halo{ii}=fitnlm(1e6*tt,yy(:,ii),mramsey_mdl,mramsey_par0,'CoefficientNames',mramsey_cname);
end

% get fit params
mramsey_fpar_halo=arrayfun(@(I) cellfun(@(f) f.Coefficients.Estimate(I),mramsey_fit_halo),...
    1:numel(mramsey_cname),'UniformOutput',false);
mramsey_fpar_halo=cat(1,mramsey_fpar_halo{:});

mramsey_fparerr_halo=arrayfun(@(I) cellfun(@(f) f.Coefficients.SE(I),mramsey_fit_halo),...
    1:numel(mramsey_cname),'UniformOutput',false);
mramsey_fparerr_halo=cat(1,mramsey_fparerr_halo{:});
% dims: PAR# X MJ

om_mramsey_halo=mramsey_fpar_halo(idx_om_mramsey,:);
omerr_mramsey_halo=mramsey_fparerr_halo(idx_om_mramsey,:);

% magnetic field
B_mramsey_halo=1e6*om_mramsey_halo/(2*pi*C_gymag);
Berr_mramsey_halo=1e6*omerr_mramsey_halo/(2*pi*C_gymag);


%%% pop-inversion
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
h=figure('Name',figname,'Units',f_units,'Position',[0.2,0.2,0.5,0.2],'Renderer',f_ren);
hold on;

pexp=ploterr(1e6*tau,P_halo_avg,...
    [],P_halo_se,'o','hhxy',0);
set(pexp(1),'MarkerFaceColor',clvir(1,:),'MarkerEdgeColor',cvir(1,:));
set(pexp(2),'Color',cvir(1,:));

tt_fit=1e6*linspace(min(tau),max(tau),1e3);     % x-axis range for fitted curve
yy_fit=feval(Pramsey_fit_halo,tt_fit);
pfit=plot(tt_fit,yy_fit,'LineStyle',line_sty{1},'Color',clvir(1,:));
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

% save fig
if do_save_figs
    savefigname=sprintf('fig_%s_%s',figname,getdatetimestr);
    fpath=fullfile(dir_save,savefigname);
    
    saveas(h,strcat(fpath,'.fig'),'fig');
    print(h,strcat(fpath,'.svg'),'-dsvg');
end


%% atom number distribution: nk_m
%%% Spatial zones
% construct spatial zones at latlon grid + solid angle
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


%% n,P distribution
%%% atom number
Nm_k=cell(numel(tau),1);     % #atoms in k,mj-zone categorise by exp param

% evaluate
for ii=1:numel(k_ramsey_phi0)
    tk=k_ramsey_phi0{ii};        % exp data for this param
    
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
%%%% spin-comp
pm_k=cellfun(@(n,N) n./N,Nm_k,N_k,'UniformOutput',false);
pm_k_avg=cellfun(@(x) squeeze(mean(x,3,'omitnan')),pm_k,'UniformOutput',false);
pm_k_std=cellfun(@(x) squeeze(std(x,0,3,'omitnan')),pm_k,'UniformOutput',false);

% multidim array: T_DELAY X AZ X EL X MJ
pm_k_avg=cell2mat(cellfun(@(a) reshape(a,[1,size(a)]),pm_k_avg,'UniformOutput',false));
pm_k_std=cell2mat(cellfun(@(a) reshape(a,[1,size(a)]),pm_k_std,'UniformOutput',false));
pm_k_se=pm_k_std./sqrt(n_shot);     


%%%% p-inversion
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
%%% spin-comp
mramsey_fit_k=cell(n_az,n_el,2);    % dim: AZ X EL X MJ
for ii=1:n_zone
    [iaz,iel]=ind2sub([n_az,n_el],ii);
    for jj=1:2
        % collate data to fit
        yy=cellfun(@(x) squeeze(x(iaz,iel,:,jj)),pm_k,'UniformOutput',false);
        yy=vertcat(yy{:});       
        
        % fit
        mramsey_fit_k{iaz,iel,jj}=fitnlm(1e6*tt,yy,mramsey_mdl,...
            mramsey_par0,'CoefficientNames',mramsey_cname);
    end
end

% get fit params
mramsey_fpar_k=arrayfun(@(I) cellfun(@(f) f.Coefficients.Estimate(I),mramsey_fit_k),...
    1:numel(mramsey_cname),'UniformOutput',false);
mramsey_fpar_k=cat(ndims(mramsey_fpar_k{1})+1,mramsey_fpar_k{:});

mramsey_fparerr_k=arrayfun(@(I) cellfun(@(f) f.Coefficients.SE(I),mramsey_fit_k),...
    1:numel(mramsey_cname),'UniformOutput',false);
mramsey_fparerr_k=cat(ndims(mramsey_fparerr_k{1})+1,mramsey_fparerr_k{:});
% dims: AZ X EL X MJ X PAR#

% get freq
om_mramsey_k=mramsey_fpar_k(:,:,:,idx_om_mramsey);
omerr_mramsey_k=mramsey_fparerr_k(:,:,:,idx_om_mramsey);
% AZ X EL X MJ

% % outlier to NaN
% b_outlier=isoutlier(om_mramsey_mode);    % find outlier by Med Abs Dev
% om_mramsey_mode_raw=om_mramsey_mode;      % store original raw data set
% omerr_mramsey_mode_raw=omerr_mramsey_mode;
% om_mramsey_mode(b_outlier)=NaN;          % outlier --> NaN
% omerr_mramsey_mode(b_outlier)=NaN;

% om_mramsey_mode_0=squeeze(mean(om_mramsey_mode,1,'omitnan'));   % mode-resolved L-freq avgd thru PHI_DELAY
% n_samp=squeeze(sum(~b_outlier,1));       % num exp samples to average over
% omerr_mramsey_mode_0=sqrt(squeeze(sum(omerr_mramsey_mode.^2,1,'omitnan')))./n_samp;

% magnetic field
B_mramsey_k=1e6*om_mramsey_k/(2*pi*C_gymag);
Berr_mramsey_k=1e6*omerr_mramsey_k/(2*pi*C_gymag);


%%% pop-inversion
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

% magnetic field
Bk_Pramsey=1e6*omk_Pramsey/(2*pi*C_gymag);
Berrk_Pramsey=1e6*omerrk_Pramsey/(2*pi*C_gymag);


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

% % outlier to NaN
% b_outlier=isoutlier(om_eq_Pramsey);    % find outlier by Med Abs Dev
% omk_Pramsey_raw=om_eq_Pramsey;      % store original raw data set
% omerr_eq_Pramsey_raw=omerr_eq_Pramsey;
% om_eq_Pramsey(b_outlier)=NaN;          % outlier --> NaN
% omerr_eq_Pramsey(b_outlier)=NaN;

% om_eq_Pramsey_0=squeeze(mean(om_eq_Pramsey,'omitnan'));   % mode-resolved L-freq avgd thru PHI_DELAY
% n_samp=squeeze(sum(~b_outlier));       % num exp samples to average over
% omerr_eq_Pramsey_0=sqrt(squeeze(sum(omerr_eq_Pramsey.^2,1,'omitnan')))./n_samp;

% magnetic field
B_Pramsey_eq=1e6*om_Pramsey_eq/(2*pi*C_gymag);
Berr_Pramsey_eq=1e6*omerr_Pramsey_eq/(2*pi*C_gymag);


%% VIS: Mode-resolved Ramsey fringe: Pm
figname='mramsey_fringe_mode';
h=figure('Name',figname,'Units',f_units,'Position',[0.2,0.2,0.5,0.2],'Renderer',f_ren);
hold on;

xx=1e6*linspace(min(tau),max(tau),1e3);     % x-axis range for fitted curve

for ii=1:naz_disp
    iaz=iaz_disp(ii);
    tp=squeeze(pm_k_avg(:,iaz,iel_0,:));   % T X MJ
    tperr=squeeze(pm_k_se(:,iaz,iel_0,:));
    
    subplot(1,naz_disp,ii);
    hold on;
    for jj=1:2
        pexp=ploterr(1e6*tau,tp(:,jj),[],tperr(:,jj),mark_typ{jj},'hhxy',0);
        set(pexp(1),'MarkerFaceColor',clvir(jj,:),'MarkerEdgeColor',cvir(jj,:),...
            'DisplayName',str_ss{jj});
        set(pexp(2),'Color',cvir(jj,:));
        
        yy=feval(mramsey_fit_k{iaz,iel_0,jj},xx);
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
    xlabel('pulse delay $\tau~(\mu s)$');
    ylabel('pop fraction $P_i$');
    xlim(1e6*[min(tau),max(tau)]+[-0.1,0.1]);
    ylim([-0.05,1.05]);
end

% save fig
if do_save_figs
    savefigname=sprintf('fig_%s_%s',figname,getdatetimestr);
    fpath=fullfile(dir_save,savefigname);
    
    saveas(h,strcat(fpath,'.fig'),'fig');
    print(h,strcat(fpath,'.svg'),'-dsvg');
end

%% VIS: Mode-resolved Ramsey fringe: P
% c_az=viridis(naz_disp);
c_az=palette(naz_disp);
cl_az=colshades(c_az);

figname='Pramsey_fringe_mode';
h=figure('Name',figname,'Units',f_units,'Position',[0.2,0.2,0.5,0.2],'Renderer',f_ren);
hold on;

xx=1e6*linspace(min(tau),max(tau),1e3);     % x-axis range for fitted curve

pleg=[];
for ii=1:naz_disp
    iaz=iaz_disp(ii);
    tp=squeeze(P_k_avg(:,iaz,iel_0));
    tperr=squeeze(P_k_se(:,iaz,iel_0));
    
    pexp=ploterr(1e6*tau,tp,[],tperr,mark_typ{ii},'hhxy',0);
    set(pexp(1),'MarkerFaceColor',cl_az(ii,:),'MarkerEdgeColor',c_az(ii,:),...
        'DisplayName',num2str(rad2deg(az(iaz))));
%         'DisplayName',num2str(az(iaz)/pi));
    set(pexp(2),'Color',c_az(ii,:));
    pleg(ii)=pexp(1);
    
    yy=feval(Pramseyk_fit{iaz,iel_0},xx);
    pfit=plot(xx,yy,'LineStyle',line_sty{ii},'Color',c_az(ii,:));
    uistack(pfit,'bottom');
    
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
end

lgd=legend(pleg);
% title(lgd,'$\theta/\pi$');
title(lgd,'$\theta$ (deg)');

% save fig
if do_save_figs
    savefigname=sprintf('fig_%s_%s',figname,getdatetimestr);
    fpath=fullfile(dir_save,savefigname);
    
    saveas(h,strcat(fpath,'.fig'),'fig');
    print(h,strcat(fpath,'.svg'),'-dsvg');
end

% %% VIS: Magnetometry: Pm
% idx_mJ=1;       % just use fit param from mJ=1
% 
% figname='halo_magnetometry_Pm';
% h=figure('Name',figname,'Units',f_units,'Position',[0.2,0.2,0.5,0.2],'Renderer',f_ren);
% 
% tp=plotFlatMapWrappedRad(gaz,gel,squeeze(B_mramsey_k(:,:,idx_mJ)),'rect','texturemap');
% 
% % annotation
% ax=gca;
% set(ax,'Layer','Top');
% box on;
% grid on;
% ax.FontSize=fontsize;
% ax.LineWidth=ax_lwidth;
% 
% axis tight;
% xlim([-180,180]);
% xticks(-180:90:180);
% yticks(-90:45:90);
% 
% xlabel('$\theta$ (deg)');
% ylabel('$\phi$ (deg)');
% 
% cbar=colorbar('eastoutside');
% cbar.TickLabelInterpreter='latex';
% cbar.Label.Interpreter='latex';
% cbar.Label.String='Magnetic field $\mathrm{B}$ (G)';
% cbar.Label.FontSize=fontsize;
% cbar.FontSize=fontsize;
% colormap('viridis');
% 
% 
% % save fig
% if do_save_figs
%     savefigname=sprintf('fig_%s_%s',figname,getdatetimestr);
%     fpath=fullfile(dir_save,savefigname);
%     
%     saveas(h,strcat(fpath,'.fig'),'fig');
%     print(h,strcat(fpath,'.png'),'-dpng','-r300');
% end


%% VIS: Magnetometry: P
figname='halo_magnetometry_P';
h=figure('Name',figname,'Units',f_units,'Position',[0.2,0.2,0.5,0.2],'Renderer',f_ren);


tp=plotFlatMapWrappedRad(gaz,gel,squeeze(Bk_Pramsey(:,:)),'rect','texturemap');
% TODO
% for rect projection, IMAGESC --> controllable x,y axis?

% annotation
ax=gca;
set(ax,'Layer','Top');
box on;
grid on;
ax.FontSize=fontsize;
ax.LineWidth=ax_lwidth;

axis tight;
xlim([-180,180]);
xticks(-180:90:180);
yticks(-90:45:90);

xlabel('$\theta$ (deg)');
ylabel('$\phi$ (deg)');

cbar=colorbar('eastoutside');
cbar.TickLabelInterpreter='latex';
cbar.Label.Interpreter='latex';
cbar.Label.String='Magnetic field $\mathrm{B}$ (G)';
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


% %% VIS: Magnetic tomography: Pm: equatorial
% figname='B_equatorial_tomography';
% h=figure('Name',figname,'Units',f_units,'Position',f_pos,'Renderer',f_ren);
% 
% Bk_eq=B_mramsey_k(:,:,iel_0,idx_mJ);      % magnetic field around equator
% Bkerr_eq=Berr_mramsey_k(:,:,iel_0,idx_mJ);      
% 
% hold on;
% pleg=NaN(npar_phidelay,1);
% for ii=1:npar_phidelay
%     tp=ploterr(rad2deg(az),Bk_eq(ii,:),[],Bkerr_eq(ii,:),'o','hhxy',0);
%     set(tp(1),'Marker',mark_typ{ii},...
%         'MarkerFaceColor',clvir2(ii,:),'MarkerEdgeColor',cvir2(ii,:),...
%         'DisplayName',num2str(rad2deg(par_comp{2}(ii)),3));
%     set(tp(2),'Color',cvir2(ii,:));
%     pleg(ii)=tp(1);
% end
% 
% % annotation
% box on;
% ax=gca;
% set(ax,'Layer','Top');
% ax.FontSize=fontsize;
% ax.LineWidth=ax_lwidth;
% 
% title('Equatorial tomography $\phi=0$');
% 
% xlabel('Azimuthal angle $\theta$ [$^\circ$]');
% ylabel('Magnetic field $\mathrm{B}$ [G]');
% 
% xlim([-180,180]);
% ax.XTick=-180:90:180;
% 
% % lgd=legend(pleg);
% 
% % save fig
% if do_save_figs
%     savefigname=sprintf('fig_%s_%s',figname,getdatetimestr);
%     fpath=fullfile(dir_save,savefigname);
%     
%     saveas(h,strcat(fpath,'.fig'),'fig');
%     print(h,strcat(fpath,'.svg'),'-dsvg');
% end


%% VIS: Magnetic tomography: equatorial: P
figname='B_equatorial_tomography_P';
% h=figure('Name',figname,'Units',f_units,'Position',f_pos,'Renderer',f_ren);
h=figure('Name',figname,'Units',f_units,'Position',[0.2,0.2,0.5,0.2],'Renderer',f_ren);
hold on;

%%% EQ-integrated 
H=shadedErrorBar([-180,180],B_Pramsey_eq*[1,1],Berr_Pramsey_eq*[1,1]);

%%% k-resolved
Bk_eq=Bk_Pramsey(:,iel_0);      % magnetic field around equator
Bkerr_eq=Berrk_Pramsey(:,iel_0);      

tp=ploterr(rad2deg(az),Bk_eq,[],Bkerr_eq,'o','hhxy',0);
set(tp(1),'Marker',mark_typ{1},...
    'MarkerFaceColor',clvir2(1,:),'MarkerEdgeColor',cvir2(1,:),...
    'DisplayName','');
set(tp(2),'Color',cvir2(1,:));

% annotation
box on;
ax=gca;
set(ax,'Layer','Top');
ax.FontSize=fontsize;
ax.LineWidth=ax_lwidth;

xlabel('$\theta$ (deg)');
ylabel('Magnetic field $\mathrm{B}$ (G)');

xlim([-180,180]);
ax.XTick=-180:90:180;

% lgd=legend(pleg);
% legend(H.mainLine,'equator integrated');

% save fig
if do_save_figs
    savefigname=sprintf('fig_%s_%s',figname,getdatetimestr);
    fpath=fullfile(dir_save,savefigname);
    
    saveas(h,strcat(fpath,'.fig'),'fig');
    print(h,strcat(fpath,'.svg'),'-dsvg');
end