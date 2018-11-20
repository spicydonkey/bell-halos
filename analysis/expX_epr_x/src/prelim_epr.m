%% Prelim multiplicative variance EPR-steering
%
% DKS
% 2018-11-19
%

%% CONFIGS
% 2018-11-19: first run X,Y,Z data (logbook 8 p93)
% fdata='C:\Users\HE BEC\Dropbox\phd\data\bell_epr_2018\proc\expX_epr_x\prelim_20181119\k_mm_20181119.mat';

% 2018-11-20: X,Z data includes other experiments (L8 p93)
fdata='C:\Users\HE BEC\Dropbox\phd\data\bell_epr_2018\proc\expX_epr_x\prelim_20181119\k_mm_comb_20181120.mat';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT
do_save_figs=true;
dir_save='C:\Users\HE BEC\Dropbox\phd\thesis\projects\maggrad+epr\epr\prelim_20181120';

% k-mode (A,B - spherically opposite)
alpha=0.08;

n_az=40;                	% equispaced bins
n_el=20;

az_disp=deg2rad([0,45,90]);     % azim sections (great circles) to display
el_disp=deg2rad([-30,0,30]);    % elev/lat zones to display

% EPR-steering
%   Note: symmetry from the post-selected spin-1/2 condition (RHS of
%   criteria)
g_inf = {-1, -1, 1};    % for Jx, Jy, Jz measurements

% vis
f_units='normalized';
f_pos=[0.2,0.2,0.2,0.3];
f_pos_wide=[0.2,0.2,0.25,0.3];
f_ren='painters';

cviridis=viridis(3);
[c,cl,cd]=palette(3);
c_gray=0.8*ones(1,3);
line_sty={'-','--',':','-.'};
mark_typ={'o','s','^','d'};
str_mm={'$x$','$y$','$z$'};
% str_ss={'$\vert\!\uparrow\rangle$','$\vert\!\downarrow\rangle$'};
% str_ss={'$m_J = 1$','$m_J = 0$'};
str_ss={'$\uparrow\uparrow$','$\downarrow\downarrow$','$\uparrow\downarrow$'};
mark_siz=7;
line_wid=1.5;
fontsize=12;
ax_lwidth=1.2;

%% load data
load(fdata);
% k_mm = {k_xx,k_yy,k_zz};

%% construct k-modes
lim_az=[0,pi];              % limit of azim angle (exc +pi)
phi_max=pi/4;
lim_el=[-phi_max,phi_max];

Vaz=linspace(lim_az(1),lim_az(2),n_az+1);
Vaz=Vaz(1:end-1);       % exclude last
Vel=linspace(lim_el(1),lim_el(2),n_el);

[vaz,vel]=ndgrid(Vaz,Vel);    % AZ-EL grid
n_zone=numel(vaz);

% special angles
[~,iaz_disp]=arrayfun(@(x) min(abs(Vaz-x)),az_disp);     % idx to displayable azim
[~,iel_disp]=arrayfun(@(x) min(abs(Vel-x)),el_disp);     % idx to displayable elev
naz_disp=length(iaz_disp);
nel_disp=length(iel_disp);
[~,iel_0]=min(abs(Vel));        % idx to ~zero elev angle (equator)


%% main
n_M=numel(g_inf);       % number of measurement configurations
% preallocate
Jm_AB=cell(n_M,n_az,n_el);
std_Jm_AB=NaN(n_M,n_az,n_el);
n_det_events=NaN(n_M,n_az,n_el);

for mm=1:n_M      % for measurement configuration
    tk=k_mm{mm};
    
    progressbar();
    for kk=1:n_zone
        [iaz,iel]=ind2sub([n_az,n_el],kk);
        
        taz=vaz(kk);
        tel=vel(kk);
        
        zone_k{1}=[taz,tel];       % az-el vector to A
        zone_k{2}=[taz+pi,-tel];   % B: sph-opposite to A
        
        % get num counts in A/B
        N=cell(1,2);        % spin-resolved numbers in A/B
        Ntot=cell(1,2);
        for ii=1:2
            [~,~,N{ii}]=cellfun(@(x) inCone(x,zone_k{ii}(1),zone_k{ii}(2),alpha),tk,'UniformOutput',false);
            N{ii}=cellfun(@(x) x,N{ii});
            
            Ntot{ii}=sum(N{ii},2);      % get total counts
        end
        
        % post-selection: single-atom detection at A/B
        b_1atom_k=cellfun(@(n) n==1, Ntot, 'UniformOutput',false);
        b_1atom=and(b_1atom_k{:});
        
        N_ps=cellfun(@(x) x(b_1atom,:),N,'UniformOutput',false);
        
        % angular momentum
        %   0.5*(N_up - N_down)
        Jm=cellfun(@(n) 0.5*(n(:,1)-n(:,2)),N_ps,'UniformOutput',false);
        
        % collective measurement
        tJm_AB=g_inf{mm}*Jm{1}+Jm{2};       % linear estimation
        tn_det_events=numel(tJm_AB);        % successfully post-selected events
        tstd_Jm_AB=std(tJm_AB);             % variation in coll meas of spin
        
        % store
        Jm_AB{mm,iaz,iel}=tJm_AB;
        n_det_events(mm,iaz,iel)=tn_det_events;
        std_Jm_AB(mm,iaz,iel)=tstd_Jm_AB;
        
        progressbar(kk/n_zone);
    end
end

%% Statistics
val_Jm_AB=unique(cat(1,Jm_AB{:}));      % unique inf coll meas outcome
% # events for each unique outcome
N_Jm_AB=arrayfun(@(j) cellfun(@(x) equalCount(x,j),Jm_AB),val_Jm_AB,'UniformOutput',false);

% summary
Navg_Jm_AB=cellfun(@(n) mean(reshape(n,n_M,[]),2),N_Jm_AB,'UniformOutput',false);
Nstd_Jm_AB=cellfun(@(n) std(reshape(n,n_M,[]),[],2),N_Jm_AB,'UniformOutput',false);
nsamp_Jm_AB=cellfun(@(n) size(reshape(n,n_M,[]),2),N_Jm_AB);

Navg_Jm_AB=cat(2,Navg_Jm_AB{:});
Nstd_Jm_AB=cat(2,Nstd_Jm_AB{:});
Nse_Jm_AB=Nstd_Jm_AB./sqrt(nsamp_Jm_AB);
% NOTE:  indexing: MM X JJ

% normalise
Ntot_Jm_AB=sum(Navg_Jm_AB,2);
navg_Jm_AB=Navg_Jm_AB./Ntot_Jm_AB;
nstd_Jm_AB=Nstd_Jm_AB./Ntot_Jm_AB;
nse_Jm_AB=Nse_Jm_AB./Ntot_Jm_AB;


%% EPR-steering
Sepr_lim=1/4;      % EPR-steering bound (spin-1/2)
Sepr=squeeze(std_Jm_AB(1,:,:).*std_Jm_AB(3,:,:));       % [Jx^B,Jz^B]
% Sepr=squeeze(std_Jm_AB(1,:,:).*std_Jm_AB(2,:,:));       % [Jx^B,Jy^B]
% Sepr=squeeze(std_Jm_AB(2,:,:).*std_Jm_AB(3,:,:));       % [Jy^B,Jz^B]

% statistics
Sepr_avg = mean(Sepr(:),'omitnan');
Sepr_std = std(Sepr(:),'omitnan');
nsamp_Sepr=sum(~isnan(Sepr(:)));
Sepr_se = Sepr_std/sqrt(nsamp_Sepr);

str_stat_epr=sprintf('%s%0.2g %s %0.1g','$\overline{\mathcal{S}}=$',...
    Sepr_avg,'$\pm$',Sepr_se);      % summary


%% VIS
%% vis: histogram of post-selected detection events
figname=sprintf('hist_numps');
h=figure('Name',figname,'Units',f_units,'Position',f_pos,'Renderer',f_ren);

hold on;
H=[];
for mm=1:3
    H(mm)=histogram(n_det_events(mm,:,:),...
        'BinWidth',1,...
        'FaceColor',cviridis(mm,:),...
        'Normalization','pdf','DisplayName',str_mm{mm});
end

% annotation
ax=gca;
box on;
set(ax,'Layer','Top');
ax.FontSize=fontsize;
% ax.LineWidth=ax_lwidth;
lgd=legend(H);
xlabel('\# events (post-selected)');
ylabel('PDF');

% save fig
if do_save_figs
    savefigname=sprintf('fig_%s_%s',figname,getdatetimestr);
    fpath=fullfile(dir_save,savefigname);
    
    saveas(h,strcat(fpath,'.fig'),'fig');
    print(h,strcat(fpath,'.svg'),'-dsvg');
end

%% vis: collective meas around different modes
figname=sprintf('coll_spin');
h=figure('Name',figname,'Units',f_units,'Position',f_pos,'Renderer',f_ren);

H=barwitherr(nse_Jm_AB',navg_Jm_AB','XData',val_Jm_AB);

for mm=1:3
    H(mm).DisplayName=str_mm{mm};
    H(mm).FaceColor=cviridis(mm,:);
end

% annotation
ax=gca;
box on;
set(ax,'Layer','Top');
ax.FontSize=fontsize;
% ax.LineWidth=ax_lwidth;
xlabel('Inferred collective spin $C_{i}$');
ylabel('PDF');
lgd=legend(H);
ax.XTick=val_Jm_AB;

% save fig
if do_save_figs
    savefigname=sprintf('fig_%s_%s',figname,getdatetimestr);
    fpath=fullfile(dir_save,savefigname);
    
    saveas(h,strcat(fpath,'.fig'),'fig');
    print(h,strcat(fpath,'.svg'),'-dsvg');
end

%% vis: histogram of inferred std
figname=sprintf('hist_infstd');
h=figure('Name',figname,'Units',f_units,'Position',f_pos,'Renderer',f_ren);

hold on;
H=[];
for mm=1:3
    H(mm)=histogram(std_Jm_AB(mm,:,:),...
        'BinWidth',0.025,...
        'FaceColor',cviridis(mm,:),...
        'Normalization','pdf','DisplayName',str_mm{mm});
end

% annotation
ax=gca;
box on;
set(ax,'Layer','Top');
ax.FontSize=fontsize;
% ax.LineWidth=ax_lwidth;
lgd=legend(H);
xlabel('Inferred uncertainty $\Delta_{\mathrm{inf}} S_{i}^{(B)}$');
ylabel('PDF');

% save fig
if do_save_figs
    savefigname=sprintf('fig_%s_%s',figname,getdatetimestr);
    fpath=fullfile(dir_save,savefigname);
    
    saveas(h,strcat(fpath,'.fig'),'fig');
    print(h,strcat(fpath,'.svg'),'-dsvg');
end

%% vis: EPR-steering
figname=sprintf('EPR');
h=figure('Name',figname,'Units',f_units,'Position',f_pos,'Renderer',f_ren);

% inversion symmetry
[vaz_sym,vel_sym,S_epr_sym]=autofill_cent_symm(vaz,vel,Sepr);

p=plotFlatMap(rad2deg(vel_sym),rad2deg(vaz_sym),S_epr_sym,'eckert4','texturemap');
colormap('viridis');

% annotation
ax=gca;
set(ax,'Layer','Top');
ax.FontSize=fontsize;
ax.LineWidth=ax_lwidth;

title(str_stat_epr);

cbar=colorbar('SouthOutside');
cbar.TickLabelInterpreter='latex';
cbar.Label.Interpreter='latex';
cbar.Label.String='$\Delta_{\mathrm{inf}} S_{z}^{(B)} \Delta_{\mathrm{inf}} S_{x}^{(B)}$';
cbar.Label.FontSize=fontsize;

% save fig
if do_save_figs
    savefigname=sprintf('fig_%s_%s',figname,getdatetimestr);
    fpath=fullfile(dir_save,savefigname);
    
    saveas(h,strcat(fpath,'.fig'),'fig');
    print(h,strcat(fpath,'.svg'),'-dsvg');
end

%% vis: histogram EPR-steering parameter
figname=sprintf('EPR_histogram');
h=figure('Name',figname,'Units',f_units,'Position',f_pos,'Renderer',f_ren);

H=histogram(Sepr(:),'BinWidth',0.025,'Normalization','pdf');
H.FaceColor=cviridis(2,:);
H.FaceAlpha=1;      % opaque
hold on;

% EPR-steering limit
ax=gca;
box on;
ax_ylim=ax.YLim;
ax_xlim=ax.XLim;
p_epr=patch([ax_xlim(1),Sepr_lim,Sepr_lim,ax_xlim(1)],...
    [ax_ylim(1),ax_ylim(1),ax_ylim(2),ax_ylim(2)],...
    c_gray,'EdgeColor','none');
uistack(p_epr,'bottom');
text(Sepr_lim,mean(ax_ylim),sprintf('EPR-steering'),...
    'HorizontalAlignment','center','VerticalAlignment','bottom',...
    'Rotation',90,'FontSize',fontsize-1);

% summary
errbar_y=mean(ax.YLim);
ebar_Sepr=errorbar(Sepr_avg,errbar_y,Sepr_std,'horizontal',...
    'LineWidth',line_wid,'Marker','o',...
    'DisplayName','$\mu \pm \sigma$');
xlim(ax_xlim);

% annotation
set(ax,'Layer','Top');
ax.FontSize=fontsize;
% ax.LineWidth=ax_lwidth;
xlabel('Inf. unc. product $\mathcal{S} = \Delta_{\mathrm{inf}} S_{z}^{(B)} \Delta_{\mathrm{inf}} S_{x}^{(B)}$');
ylabel('PDF');
% lgd=legend(ebar_Sepr);    
ylim(ax_ylim);

% save fig
if do_save_figs
    savefigname=sprintf('fig_%s_%s',figname,getdatetimestr);
    fpath=fullfile(dir_save,savefigname);
    
    saveas(h,strcat(fpath,'.fig'),'fig');
    print(h,strcat(fpath,'.svg'),'-dsvg');
end


%%
% %% vis
% %% observer
% h=figure;
%
% subplot(1,2,1);
% hold on;
% H(1)=histogram(n_A(:,1),'DisplayName','$\uparrow$');
% H(2)=histogram(n_A(:,2),'DisplayName','$\downarrow$');
%
% xlabel('$N_{m}$');
% ylabel('\# observations');
% lgd=legend(H);
% title('Alice');
%
% subplot(1,2,2);
% hold on;
% H(1)=histogram(n_B(:,1),'DisplayName','$\uparrow$');
% H(2)=histogram(n_B(:,2),'DisplayName','$\downarrow$');
%
% xlabel('$N_{m}$');
% ylabel('\# observations');
% lgd=legend(H);
% title('Bob');
%
% %% angular momentum
% h=figure;
% hold on;
%
% H(1)=histogram(Jn_A,'DisplayName','Alice');
% H(2)=histogram(Jn_B,'DisplayName','Bob');
%
% xlabel('Angular momentum $J_{i}$');
% ylabel('\# observations');
% lgd=legend(H);
%
% %% tot atoms
% h=figure;
% hold on;
%
% H(1)=histogram(ntot_A,'DisplayName','Alice');
% H(2)=histogram(ntot_B,'DisplayName','Bob');
%
% xlabel('$N_{\uparrow} + N_{\downarrow}$');
% ylabel('\# observations');
% lgd=legend(H);
