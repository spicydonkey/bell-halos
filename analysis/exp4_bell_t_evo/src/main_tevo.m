%% Time evolution of Psi+ halo
% DKS
% 2018-10-25

tic

%% CONFIGS
% data file
fdata='C:\Users\HE BEC\Dropbox\PhD\projects\halo_metrology\analysis\exp4_tevo\exp4_20181029.mat';

% capture region - "central symmetry"!
alpha=pi/8;        % cone half-angle
lim_az=[0,pi];              % limit of azim angle (exc +pi)
phi_max=pi/4;
lim_el=[-phi_max,phi_max];

n_az=10;                	% equispaced bins
n_el=5;

az_disp=deg2rad([0,45,90,135]);     % azim sections (great circles) to display

% g2
n_dk=7;     %15
lim_dk=[-0.2,0.2];
% n_dk=7;
% lim_dk=[-0.1,0.1];

% bootstrapping
bs_frac=1;
bs_nrep=50;

% He*
mu_B=1.3996e6;      % bohr magneton [Hz/G]
g_He=2;             % Lande g-factor He*
C_gymag=g_He*mu_B;     % He* gyromagnetic ratio (gamma) [Hz/G]

% vis
f_units='normalized';
f_pos=[0.2,0.2,0.2,0.3];
f_pos_wide=[0.2,0.2,0.25,0.3];
f_ren='painters';

[c,cl,cd]=palette(3);
% c_gray=0.6*ones(1,3);
line_sty={'-','--',':','-.'};
mark_typ={'o','s','^','d'};
% str_ss={'$\vert\!\uparrow\rangle$','$\vert\!\downarrow\rangle$'};
% str_ss={'$m_J = 1$','$m_J = 0$'};
% str_ss={'$\uparrow\uparrow$','$\downarrow\downarrow$','$\uparrow\downarrow$'};
mark_siz=7;
line_wid=1.5;
fontsize=12;
ax_lwidth=1.2;

%% load data
load(fdata);

%% PROTOTYPE: reduce data
warning('Reducing data for speed.');
k_tau{1}=k_tau{1}(1:3000,:);

% warning('DEBUG ONLY!!!!');
% k_tau=cellfun(@(k) k(1:100,:),k_tau,'UniformOutput',false);


%% create halo sections to investigate
Vaz=linspace(lim_az(1),lim_az(2),n_az+1);
Vaz=Vaz(1:end-1);       % exclude last
Vel=linspace(lim_el(1),lim_el(2),n_el);

[vaz,vel]=ndgrid(Vaz,Vel);    % AZ-EL grid
n_zone=numel(vaz);

% special angles
[~,iaz_disp]=arrayfun(@(x) min(abs(Vaz-x)),az_disp);     % idx to ~displayable azim angles
[~,iel_0]=min(abs(Vel));        % idx to ~zero elev angle (equator)

%% create g2 params
dk_ed_vec=linspace(lim_dk(1),lim_dk(2),n_dk+1);
dk_cent_vec=dk_ed_vec(1:end-1)+0.5*diff(dk_ed_vec);
[~,idx_dk0]=min(abs(dk_cent_vec));    % idx bin nearest zero
dk_ed={dk_ed_vec,dk_ed_vec,dk_ed_vec};
dk_grid_size=cellfun(@(k) length(k)-1,dk_ed);

%% MAIN
n_tau=length(k_tau);

% preallocate
c4=cell(3,n_tau,n_az,n_el);     % DIM: corr X tau X theta X phi
c3=cell(n_tau,n_az,n_el);       % DIM: tau X theta X phi
m4=NaN(3,n_tau,n_az,n_el);    
m3=NaN(n_tau,n_az,n_el);

g2=c4;
g0=m4;    
B=m3;
B0=m3;

g2_bs_mu=c4;
g2_bs_se=c4;
g0_bs_mu=m4;
g0_bs_se=m4;
B_bs_mu=m3;
B0_bs_mu=m3;
B_bs_se=m3;
B0_bs_se=m3;

for ii=1:n_tau
    %%
    k=k_tau{ii};
    n_data_size=size(k,1);
    
    progressbar(0);
    for jj=1:n_zone
        %%
        [iaz,iel]=ind2sub([n_az,n_el],jj);
        taz=vaz(jj);
        tel=vel(jj);
        
        % capture vecs in section
%         [~,b_A]=cellfun(@(x) inCone(x,taz,tel,alpha),k,'UniformOutput',false);     % in k
%         [~,b_B]=cellfun(@(x) inCone(x,taz+pi,-tel,alpha),k,'UniformOutput',false);    % in -k
%         b_AB=cellfun(@(b1,b2) b1|b2,b_A,b_B,'UniformOutput',false);        % in both
%         k_in=cellfun(@(x,b) x(b,:),k,b_AB,'UniformOutput',false);
        k_in=cellfun(@(x) inDoubleCone(x,taz,tel,alpha),k,'UniformOutput',false);

        %% g2
        tg2=summary_disthalo_g2(k_in,dk_ed,0,0,0,0);
        
        tg0=cellfun(@(g) g(idx_dk0,idx_dk0,idx_dk0),tg2);    % BB corr strength
        [tB,tB0]=g2toE(mean(tg0(1:2)),tg0(3));      % correlator
        
        %%% store
        g2(:,ii,iaz,iel)=tg2(:);
        g0(:,ii,iaz,iel)=tg0;
        B(ii,iaz,iel)=tB;
        B0(ii,iaz,iel)=tB0;
        
        %% BOOTSTRAPPING
        % set up
        bs_nsamp=round(bs_frac*n_data_size);    % # data to sample for bs
        bs_Isamp=cellfun(@(c) randi(n_data_size,[bs_nsamp,1]), cell(bs_nrep,1),...
            'UniformOutput',false);

        %%% run
        % g2
        tg2_bs=cellfun(@(I) summary_disthalo_g2(k_in(I,:),dk_ed,0,0,0,0),bs_Isamp,...
            'UniformOutput',false);
        tg2_bs=cat(1,tg2_bs{:});    % g2 dist from bs
        tg0_bs=cellfun(@(g) g(idx_dk0,idx_dk0,idx_dk0),tg2_bs);    % corr strength
        
        % B
        [tB_bs,tB0_bs]=g2toE(mean(tg0_bs(:,1:2),2),tg0_bs(:,3));
        
        %%% statistics
        tg2_bs_cat=cell(1,3);       % tidy bs sampled g2
        for kk=1:3
            tg2_bs_cat{kk}=cat(4,tg2_bs{:,kk});
        end
        tg2_bs_mu=cellfun(@(x) mean(x,4),tg2_bs_cat,'UniformOutput',false);
        tg2_bs_se=cellfun(@(x) sqrt(bs_frac)*std(x,0,4),tg2_bs_cat,'UniformOutput',false);

        tg0_bs_mu=mean(tg0_bs,1);
        tg0_bs_se=sqrt(bs_frac)*std(tg0_bs,0,1);
        
        tB_bs_mu=mean(tB_bs,1);
        tB0_bs_mu=mean(tB0_bs,1);
        tB_bs_se=sqrt(bs_frac)*std(tB_bs,0,1);
        tB0_bs_se=sqrt(bs_frac)*std(tB0_bs,0,1);
        
        %%% store
        g2_bs_mu(:,ii,iaz,iel)=tg2_bs_mu(:);
        g2_bs_se(:,ii,iaz,iel)=tg2_bs_se(:);
        g0_bs_mu(:,ii,iaz,iel)=tg0_bs_mu;
        g0_bs_se(:,ii,iaz,iel)=tg0_bs_se;
        B_bs_mu(ii,iaz,iel)=tB_bs_mu;
        B0_bs_mu(ii,iaz,iel)=tB0_bs_mu;
        B_bs_se(ii,iaz,iel)=tB_bs_se;
        B0_bs_se(ii,iaz,iel)=tB0_bs_se;
        
        progressbar(jj/n_zone);
    end
end

%% Fit time-evolution model
% simplest model: time-independent asymmetry of Bell triplet
% B_pi/2 = cos(hbar*gamma* deltaB * t)

% set up model and solver
mdl_tevo.mdl='y~cos(om*x)';
mdl_tevo.cname={'om'};
mdl_tevo.fopt=statset('TolFun',10^-10,'TolX',10^-10,'MaxIter',10^6,'UseParallel',0);
par0=1;

% preproc
T_so=0.8;           % MODEL PARAM: asymmetry "switch-on" time = 0.8 ms
tau0=tau-T_so;      % time evolution since switch-on
% tau0_fit=linspace(min(tau0),max(tau0),1e3);     % tau vals to eval fit
tau0_fit=linspace(0,2,1e3);     % tau vals to eval fit

% fit
mdl_tevo.fit=cell(n_az,n_el);
for ii=1:n_zone
    [iaz,iel]=ind2sub([n_az,n_el],ii);
    
    tB0=B0(:,iaz,iel);
    mdl_tevo.fit{iaz,iel}=fitnlm(tau0,tB0,mdl_tevo.mdl,par0,...
        'CoefficientNames',mdl_tevo.cname,'Options',mdl_tevo.fopt);
end

% eval fitted model
mdl_tevo.fit_par=cellfun(@(m) m.Coefficients.Estimate,mdl_tevo.fit);
mdl_tevo.fit_par_se=cellfun(@(m) m.Coefficients.SE,mdl_tevo.fit);
B0_fit=cellfun(@(f) feval(f,tau0_fit),mdl_tevo.fit,'UniformOutput',false);

% diff B field 
om_fit=1e3*mdl_tevo.fit_par;        % fitted omega [rad/s]
om_se_fit=1e3*mdl_tevo.fit_par_se;
deltaB=om_fit/(2*pi*C_gymag);        % diff in B-field strength [G]
deltaB_se=om_se_fit/(2*pi*C_gymag);   % standard error (fit)

%% vis: corrected parity - ALL
%%% ALL
[cc,ccl,ccd]=palette(n_zone);

h=figure('Name','triplet_halo_tevo','Units',f_units,'Position',f_pos,'Renderer',f_ren);
hold on;

for ii=1:n_zone
    [iaz,iel]=ind2sub([n_az,n_el],ii);
    
    tp=ploterr(tau,B0(:,iaz,iel),[],B0_bs_se(:,iaz,iel),'-o','hhxy',0);
    set(tp(1),'MarkerSize',mark_siz,'LineWidth',line_wid,...
        'MarkerFaceColor',ccl(ii,:),'Color',cc(ii,:));
    set(tp(2),'LineWidth',line_wid,'Color',cc(ii,:));
    
end

% annotate
box on;
ax=gca;
set(ax,'Layer','Top');
xlabel('$\tau~[\textrm{ms}]$');
ylabel('Parity $\bar{\mathcal{B}}_{\pi/2}$');
ax.FontSize=fontsize;
ax.LineWidth=ax_lwidth;
ylim([-1.2,1.2]);

%% vis: Polar distribution
for ii=1:length(iaz_disp)
    iaz=iaz_disp(ii);       % azim idx to displayable great circle
    taz=Vaz(iaz);
    figname=sprintf('B0_tevo_polar_%0.0f',rad2deg(taz));
    
    % figure
    figure('Name',figname,...
        'Units',f_units,'Position',f_pos_wide,'Renderer',f_ren);
    hold on;
    
    pleg=NaN(n_el,1);
    for jj=1:n_el
        tp=ploterr(tau,squeeze(B0(:,iaz,jj)),[],squeeze(B0_bs_se(:,iaz,jj)),'o','hhxy',0);
        set(tp(1),'MarkerSize',mark_siz,'LineWidth',line_wid,...
            'MarkerFaceColor',ccl(jj,:),'Color',cc(jj,:),'DisplayName',num2str(rad2deg(Vel(jj)),2));
        set(tp(2),'LineWidth',line_wid,'Color',cc(jj,:));
        pleg(jj)=tp(1);
        
        % fitted model
        tpf=plot(tau0_fit+T_so,B0_fit{iaz,jj},'-',...
            'LineWidth',line_wid,'Color',cc(jj,:));
        uistack(tpf,'bottom');
    end
    
    % annotation
    box on;
    ax=gca;
    set(ax,'Layer','Top');
    xlabel('$\tau~[\textrm{ms}]$');
    ylabel('Parity $\bar{\mathcal{B}}_{\pi/2}$');
    ax.FontSize=fontsize;
    ax.LineWidth=1.2;
    ylim([-1.2,1.2]);
    titlestr=sprintf('%s %0.0f','$\theta=$',rad2deg(taz));
    title(titlestr);
    lgd=legend(pleg,'Location','EastOutside');
    title(lgd,'Latitude $\phi$ (deg)');
end

%% vis: Equatorial distribution
figname='B0_tevo_eqt';
figure('Name',figname,...
    'Units',f_units,'Position',[0.2,0.2,0.5,0.2],'Renderer',f_ren);
hold on;

disp_iaz=iaz_disp;

% pleg=NaN(n_az,1);
pleg=NaN(numel(disp_iaz),1);
% for ii=1:n_az
for ii=1:numel(disp_iaz)
    tiaz=disp_iaz(ii);
    tp=ploterr(tau,squeeze(B0(:,tiaz,iel_0)),[],squeeze(B0_bs_se(:,tiaz,iel_0)),'o','hhxy',0);
%     tp=ploterr(tau,squeeze(B0(:,ii,iel_0)),[],squeeze(B0_bs_se(:,ii,iel_0)),'o','hhxy',0);
    set(tp(1),'Marker',mark_typ{ii},'MarkerSize',mark_siz,'LineWidth',line_wid,...
        'MarkerFaceColor',ccl(ii,:),'Color',cc(ii,:),'DisplayName',num2str(rad2deg(Vaz(tiaz)),2));
    set(tp(2),'LineWidth',line_wid,'Color',cc(ii,:));
    pleg(ii)=tp(1);
    
    %%% fitted model
    % stationary before tau0 (Psi+ stationary --> parity=1)
    tpf=plot([0,T_so],[1,1],...
        line_sty{ii},'LineWidth',line_wid,'Color',cc(ii,:));
    uistack(tpf,'bottom');
    
    % dynamics after turn-on
    tpf=plot(tau0_fit+T_so,B0_fit{tiaz,iel_0},...
        line_sty{ii},'LineWidth',line_wid,'Color',cc(ii,:));
    uistack(tpf,'bottom');
end

% annotation
box on;
ax=gca;
set(ax,'Layer','Top');
xlabel('$\tau~(\textrm{ms})$');
% ylabel('Parity $\bar{\mathcal{B}}_{\pi/2}$');
ylabel('parity');
ax.FontSize=fontsize;
ax.LineWidth=1.2;
xlim([0.7,1.8]);
ylim([-1.2,1.2]);

% lgd=legend(pleg,'Location','SouthWest');
% title(lgd,'Azimuth $\theta$ (deg)');


%% vis: B0 distribution (3D)
H=[];
for ii=1:n_tau
    %%
    figname=sprintf('B0_sphdist_3d_%0.2f',tau(ii));
    H(ii)=figure('Name',figname,'Units',f_units,'Position',f_pos,'Renderer',f_ren);
    tB0=squeeze(B0(ii,:,:));

%     tp=plot_sph_surf(vaz,vel,tB0);
    [vazf,velf,tB0f]=autofill_cent_symm(vaz,vel,tB0);
    tp=plot_sph_surf(vazf,velf,tB0f);    

    % annotation
%     tp.FaceAlpha=1;

    titlestr=sprintf('%s %0.2f ms','$\tau=$',tau(ii));
    title(titlestr);

    ax=gca;
    ax.FontSize=fontsize;
    
    axis on;
    box on;
    xlabel('$k_x$');
    ylabel('$k_y$');
    zlabel('$k_z$');
    xlim([-1,1]);
    ylim([-1,1]);
    zlim([-1,1]);
    
    cbar=colorbar();
    caxis([-1,1]);
    cbar.Limits=[-1,1];
    cbar.TickLabelInterpreter='latex';
    cbar.Label.Interpreter='latex';
    cbar.Label.String='Parity $\bar{\mathcal{B}}_{\pi/2}$';
    cbar.Label.FontSize=fontsize;
end

%% vis: B0 distribution (2D)
H=[];
for ii=1:n_tau
    %%
    figname=sprintf('B0_sphdist_2d_%0.2f',tau(ii));
    H(ii)=figure('Name',figname,'Units',f_units,'Position',f_pos,'Renderer',f_ren);
    tB0=squeeze(B0(ii,:,:));
    
%     tp=plotFlatMap(rad2deg(vel),rad2deg(vaz),tB0,'eckert4','texturemap');
    [vazf,velf,tB0f]=autofill_cent_symm(vaz,vel,tB0);
    tp=plotFlatMap(rad2deg(velf),rad2deg(vazf),tB0f,'eckert4','texturemap');
    
    % annotation
    titlestr=sprintf('%s %0.2f ms','$\tau=$',tau(ii));
    title(titlestr);

    ax=gca;
    ax.FontSize=fontsize;
    
    cbar=colorbar('southoutside');
    caxis([-1,1]);
    cbar.Limits=[-1,1];
    cbar.TickLabelInterpreter='latex';
    cbar.Label.Interpreter='latex';
    cbar.Label.String='Parity $\bar{\mathcal{B}}_{\pi/2}$';
    cbar.Label.FontSize=fontsize;
end

%% vis: deltaB (asymmetry measure) around halo: t-indep model (3D SPH)
h=figure('Name','deltaB_sphdist_3d','Units',f_units,'Position',f_pos,'Renderer',f_ren);

% tp=plot_sph_surf(vaz,vel,1e3*deltaB);
[vazf,velf,deltaBf]=autofill_cent_symm(vaz,vel,deltaB);
tp=plot_sph_surf(vazf,velf,1e3*deltaBf);

% annotation
ax=gca;
ax.FontSize=fontsize;

axis on;
box on;
xlabel('$k_x$');
ylabel('$k_y$');
zlabel('$k_z$');
xlim([-1,1]);
ylim([-1,1]);
zlim([-1,1]);

cbar=colorbar();
cbar.TickLabelInterpreter='latex';
cbar.Label.Interpreter='latex';
cbar.Label.String='$\Delta \mathrm{B}$ [mG]';
cbar.Label.FontSize=fontsize;

%% vis: deltaB (asymmetry measure) around halo: t-indep model (2D PROJ MAP)
h=figure('Name','deltaB_sphdist_2d','Units',f_units,'Position',f_pos,'Renderer',f_ren);

% tp=plotFlatMap(rad2deg(vel),rad2deg(vaz),1e3*deltaB,'rect','texturemap');
[vazf,velf,deltaBf]=autofill_cent_symm(vaz,vel,deltaB);
tp=plotFlatMap(rad2deg(velf),rad2deg(vazf),1e3*deltaBf,'eckert4','texturemap');

% annotation
ax=gca;
ax.FontSize=fontsize;

cbar=colorbar('SouthOutside');
cbar.TickLabelInterpreter='latex';
cbar.Label.Interpreter='latex';
cbar.Label.String='$\Delta \mathrm{B}$ [mG]';
cbar.Label.FontSize=fontsize;

%% vis: deltaB (rect)
h=figure('Name','deltaB_sphdist_2d','Units',f_units,'Position',[0.2,0.2,0.5,0.2],'Renderer',f_ren);

[vazf,velf,deltaBf]=autofill_cent_symm(vaz,vel,deltaB);
tp=plotFlatMapWrappedRad(vazf,velf,1e3*deltaBf,'rect','texturemap');

% label ROI
hold on;
for ii=1:numel(disp_iaz)
    tiaz=disp_iaz(ii);
    tp=plot(rad2deg(Vaz(tiaz)),rad2deg(Vel(iel_0)),'Marker',mark_typ{ii},...
        'MarkerEdgeColor',cc(ii,:),'MarkerFaceColor',ccl(ii,:),...
        'MarkerSize',10,'LineWidth',line_wid);
end

% annotation
ax=gca;
set(ax,'Layer','Top');
box on;
grid on;
ax.FontSize=fontsize;
ax.LineWidth=ax_lwidth;

axis tight;
% xlim([-180,180]);
xticks(-180:90:180);
yticks(-90:45:90);

xlabel('$\theta$ (deg)');
ylabel('$\phi$ (deg)');

cbar=colorbar('eastoutside');
cbar.TickLabelInterpreter='latex';
cbar.Label.Interpreter='latex';
cbar.Label.String='$\Delta \mathrm{B}$ (mG)';
cbar.Label.FontSize=fontsize;
cbar.FontSize=fontsize;


%% VIS: Equatorial tomography
%%% configs
p_xlim=[0,180];     % periodic BC

ctheme=magma(5);
ctheme=ctheme(1:end-1,:);
cltheme=colshades(ctheme);

idx_col=4;      

%%% DATA
dB_eq=deltaB(:,iel_0);      % delta magnetic field around equator
dBerr_eq=deltaB_se(:,iel_0);      

% fit
sine_mdl='y~c+amp*cos(2*x+phi)';      % periodic boundary condition fixes OMEGA
sine_cname={'amp','phi','c'};
sine_par0=[0.05,0,0.1];

fit_dB_eq=fitnlm(Vaz,1e3*dB_eq,sine_mdl,sine_par0,'CoefficientNames',sine_cname);
xx=linspace(0,pi,1e3);
yy=feval(fit_dB_eq,xx);

%%% PLOT
figname='dB_equatorial_tomography';
h=figure('Name',figname,'Units',f_units,'Position',[0.2,0.2,0.5,0.2],'Renderer',f_ren);

hold on;

% data
tp=ploterr(rad2deg(Vaz),1e3*dB_eq,[],1e3*dBerr_eq,'o','hhxy',0);
set(tp(1),'Marker','o','LineWidth',line_wid,...
    'MarkerFaceColor',cltheme(idx_col,:),'MarkerEdgeColor',ctheme(idx_col,:),...
    'DisplayName','');
set(tp(2),'Color',ctheme(idx_col,:),'LineWidth',line_wid);
pleg=tp(1);

% fit
pfit=plot(rad2deg(xx),yy,'LineStyle','-','Color',ctheme(idx_col,:),'LineWidth',line_wid);
uistack(pfit,'bottom');

% annotation
box on;
ax=gca;
set(ax,'Layer','Top');
ax.FontSize=fontsize;
ax.LineWidth=ax_lwidth;

xlabel('Azimuthal angle $\theta$ (deg)');
ylabel('$\Delta \mathrm{B}$ (mG)');

ylim_0=ax.YLim;
ylim([0,ylim_0(2)]);
xlim(p_xlim);
ax.XTick=0:45:180;

% lgd=legend(pleg);

%% End of script
toc