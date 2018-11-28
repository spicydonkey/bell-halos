%% ANALYSIS: integrated g2(dk) (cartesian) and "uniform" alignment
% DKS
% 2018-11-06

tic

%% CONFIGS
% data file
%   contains spin,momentum resolved vectors
fdata='C:\Users\HE BEC\Dropbox\phd\data\bell_epr_2018\proc\exp4_tevo\exp4_20181029.mat';

% save
do_save_figs=true;
dir_save='C:\Users\HE BEC\Dropbox\PhD\thesis\projects\maggrad\Gint_unif_20181114';

% k-modes
alpha=pi/8;                 % cone half-angle
lim_az=[0,pi];              % limit of azim angle (exc +pi)
phi_max=pi/4;
lim_el=[-phi_max,phi_max];
n_az=24;                	% equispaced bins
n_el=15;
az_disp=deg2rad([0,45,90]);     % azim sections (great circles) to display
el_disp=deg2rad([-30,0,30]);    % elev/lat zones to display

% g2 
n_dk=15;
lim_dk=[-0.2,0.2];

% integration
dk_int_lim=0.1;

% bootstrapping
bs_frac=0.5;
bs_nrep=10;

% FIT
do_fit=true;

% He*
C_gymag=2.8e6;     % gyromagnetic ratio (gamma) [Hz/G]

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
str_ss={'$\uparrow\uparrow$','$\downarrow\downarrow$','$\uparrow\downarrow$'};
mark_siz=7;
line_wid=1.5;
fontsize=12;
ax_lwidth=1.2;

%% configure g2
dk_ed_vec=linspace(lim_dk(1),lim_dk(2),n_dk+1);
dk_cent_vec=dk_ed_vec(1:end-1)+0.5*diff(dk_ed_vec);
[~,idx_dk0]=min(abs(dk_cent_vec));    % idx bin nearest zero
dk_ed={dk_ed_vec,dk_ed_vec,dk_ed_vec};
dk_grid_size=cellfun(@(k) length(k)-1,dk_ed);

% integrable region config
dIk_int=round(dk_int_lim/(range(lim_dk)/n_dk));
intFilter=@(x) x(idx_dk0-dIk_int:idx_dk0+dIk_int,...
    idx_dk0-dIk_int:idx_dk0+dIk_int,idx_dk0-dIk_int:idx_dk0+dIk_int);

%% load data
load(fdata);

%% PROTOTYPE: reduce data
warning('Reducing data for speed.');
k_tau{1}=k_tau{1}(1:3000,:);

% warning('DEBUG ONLY!!!!');
% k_tau=cellfun(@(k) k(1:100,:),k_tau,'UniformOutput',false);

% % TEST K-ALIGNMENT SENSITIVITY
% IDX_TAU=4;
% k_tau={k_tau{IDX_TAU}};
% tau=tau(IDX_TAU);

%% k-modes
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
n_tau=length(k_tau);

G_int=NaN(n_tau,n_az,n_el,3);
G_sum=NaN(n_tau,n_az,n_el);
B=NaN(n_tau,n_az,n_el);
B0=NaN(n_tau,n_az,n_el);

G_int_se=NaN(n_tau,n_az,n_el,3);
G_sum_se=NaN(n_tau,n_az,n_el);
B_se=NaN(n_tau,n_az,n_el);
B0_se=NaN(n_tau,n_az,n_el);

for ii=1:n_tau
    k=k_tau{ii};
    n_data_size=size(k,1);
    
    progressbar(0);
    for jj=1:n_zone
        %%
        [iaz,iel]=ind2sub([n_az,n_el],jj);
        taz=vaz(jj);
        tel=vel(jj);
        
        k_cone=cellfun(@(x) inDoubleCone(x,taz,tel,alpha),k,'UniformOutput',false);
        
        %% g2
        tg2=summary_disthalo_g2(k_cone,dk_ed,0,0,0,0);
        
        % volume integrated G2
        tG_int=cellfun(@(g) meanall(intFilter(g),'omitnan'),tg2);
        tG_sum=0.5*sum(tG_int([1 2 3 3])-1);
        [tB,tB0]=g_to_E(tG_int);
    
        %% BOOTSTRAPPING
        % set up
        bs_nsamp=round(bs_frac*n_data_size);    % # data to sample for bs
        bs_Isamp=cellfun(@(c) randi(n_data_size,[bs_nsamp,1]), cell(bs_nrep,1),...
            'UniformOutput',false);

        %%% run
        % g2
        tg2_bs=cellfun(@(I) summary_disthalo_g2(k_cone(I,:),dk_ed,0,0,0,0),bs_Isamp,...
            'UniformOutput',false);
        tg2_bs=cat(1,tg2_bs{:});    % g2 dist from bs
        
        % volume integrated G2
        tG_int_bs=cellfun(@(g) meanall(intFilter(g),'omitnan'),tg2_bs);
        tG_sum_bs=0.5*sum([tG_int_bs,tG_int_bs(:,3)],2)-2;
        [tB_bs,tB0_bs]=g_to_E(tG_int_bs);
        
        % stat
        tG_int_mu=mean(tG_int_bs,1);
        tG_int_se=sqrt(bs_frac)*std(tG_int_bs,[],1);
        tG_sum_mu=mean(tG_sum_bs,1);
        tG_sum_se=sqrt(bs_frac)*std(tG_sum_bs,[],1);
        tB_mu=mean(tB_bs,1);
        tB_se=sqrt(bs_frac)*std(tB_bs,[],1);
        tB0_mu=mean(tB0_bs,1);
        tB0_se=sqrt(bs_frac)*std(tB0_bs,[],1);
                
        %% store
        G_int(ii,iaz,iel,:)=tG_int;
        G_sum(ii,iaz,iel)=tG_sum;
        B(ii,iaz,iel)=tB;
        B0(ii,iaz,iel)=tB0;
        
        G_int_se(ii,iaz,iel,:)=tG_int_se;
        G_sum_se(ii,iaz,iel)=tG_sum_se;
        B_se(ii,iaz,iel)=tB_se;
        B0_se(ii,iaz,iel)=tB0_se;
        
        progressbar(jj/n_zone);
    end
end

%% ANALYSIS: Magnetic gradiometry
if do_fit
    %% fit time-evolution model
    % simplest model: time-independent asymmetry of Bell triplet
    % B_pi/2 = cos(hbar*gamma* deltaB * t)

    % set up model and solver
    mdl_tevo.mdl='y~cos(om*x)';
    mdl_tevo.cname={'om'};
    mdl_tevo.fopt=statset('TolFun',10^-10,'TolX',10^-10,'MaxIter',10^6,'UseParallel',0);
    par0=1;
    
    % preproc
    tau0=tau-tau(1);        % time evolution since T(PSI+)~0.8 [ms]
    tau0_fit=linspace(min(tau0),max(tau0),1e3);
    
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
    
    %% cull outliers
    % get outliers: fit uncertainty more than 1sig from median    
    % A) # stdevs
%     n_sig_outlier=0.3;
%     dx_outlier=n_sig_outlier*std(deltaB_se(:));
    % B) hard limit
%     dx_outlier=1e-3/(2*pi)^2;        % hard bound
    dx_outlier=Inf;
    
%     b_fit_outlier=abs(deltaB_se-median(deltaB_se(:)))>dx_outlier;    
    b_fit_outlier=isoutlier(deltaB_se);

    % outliers to NaN
    deltaB_filt=deltaB;
    deltaB_filt(b_fit_outlier)=NaN;
    
    deltaB_se_filt=deltaB_se;
    deltaB_se_filt(b_fit_outlier)=NaN;
end

%% DATA VISUALISATION
[cc,ccl,ccd]=palette(n_zone);

%% vis: Polar distribution
for ii=1:naz_disp
    iaz=iaz_disp(ii);       % azim idx to displayable great circle
    taz=Vaz(iaz);
    figname=sprintf('B0_tevo_polar_%0.0f',rad2deg(taz));
    
    % figure
    h=figure('Name',figname,...
        'Units',f_units,'Position',f_pos_wide,'Renderer',f_ren);
    hold on;
    
    pleg=NaN(nel_disp,1);
    for jj=1:nel_disp
        iel=iel_disp(jj);
        
        tp=ploterr(tau,squeeze(B0(:,iaz,iel)),[],squeeze(B0_se(:,iaz,iel)),'o','hhxy',0);
        set(tp(1),'MarkerSize',mark_siz,'LineWidth',line_wid,...
            'MarkerFaceColor',ccl(jj,:),'Color',cc(jj,:),'DisplayName',num2str(rad2deg(Vel(iel)),2));
        set(tp(2),'LineWidth',line_wid,'Color',cc(jj,:));
        pleg(jj)=tp(1);
        
        if do_fit
            % fitted model
            tpf=plot(tau0_fit+tau(1),B0_fit{iaz,iel},'-',...
                'LineWidth',line_wid,'Color',cc(jj,:));
            uistack(tpf,'bottom');
        end
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
    
    % save fig
    if do_save_figs
        savefigname=sprintf('fig_%s_%s',h.Name,getdatetimestr);
        fpath=fullfile(dir_save,savefigname);
        
        saveas(h,strcat(fpath,'.fig'),'fig');
        print(h,strcat(fpath,'.svg'),'-dsvg');
    end
end

%% vis: Equatorial distribution
figname='B0_tevo_eqt';
h=figure('Name',figname,...
    'Units',f_units,'Position',f_pos_wide,'Renderer',f_ren);
hold on;

pleg=NaN(naz_disp,1);
for ii=1:naz_disp
    iaz=iaz_disp(ii);
    tp=ploterr(tau,squeeze(B0(:,iaz,iel_0)),[],squeeze(B0_se(:,iaz,iel_0)),'o','hhxy',0);
    set(tp(1),'MarkerSize',mark_siz,'LineWidth',line_wid,...
        'Marker',mark_typ{mod(ii,length(mark_typ))+1},...     % will be shifted by 1
        'MarkerFaceColor',ccl(ii,:),'Color',cc(ii,:),'DisplayName',num2str(rad2deg(Vaz(iaz)),2));
    set(tp(2),'LineWidth',line_wid,'Color',cc(ii,:));
    pleg(ii)=tp(1);

    if do_fit
        % fitted model
        tpf=plot(tau0_fit+tau(1),B0_fit{iaz,iel_0},'-',...
            'LineWidth',line_wid,'Color',cc(ii,:));
        uistack(tpf,'bottom');
    end
end

% annotation
title('Equatorial');
box on;
ax=gca;
set(ax,'Layer','Top');
xlabel('$\tau~[\textrm{ms}]$');
ylabel('Parity $\bar{\mathcal{B}}_{\pi/2}$');
ax.FontSize=fontsize;
ax.LineWidth=1.2;
ylim([-1.2,1.2]);
titlestr=sprintf('Equator (%s %0.2g)','$\phi=$',Vel(iel_0));
title(titlestr);
lgd=legend(pleg,'Location','EastOutside');
title(lgd,'Azimuth $\theta$ (deg)');

% save fig
if do_save_figs
    savefigname=sprintf('fig_%s_%s',h.Name,getdatetimestr);
    fpath=fullfile(dir_save,savefigname);
    
    saveas(h,strcat(fpath,'.fig'),'fig');
    print(h,strcat(fpath,'.svg'),'-dsvg');
end

%% deltaB distribution
%% MEAN: deltaB (asymmetry measure) around halo: t-indep model (2D PROJ MAP)
if do_fit
    h=figure('Name','deltaB_sphdist_2d','Units',f_units,'Position',f_pos,'Renderer',f_ren);
    
    [vazf,velf,deltaBf]=autofill_cent_symm(vaz,vel,deltaB_filt);
    tp=plotFlatMap(rad2deg(velf),rad2deg(vazf),1e3*deltaBf,'eckert4','texturemap');
    
    % % rectangle
    % tp=plotFlatMap(rad2deg(velf),rad2deg(vazf),1e3*deltaBf,'rect','texturemap');
    % xlim([-180,180]);
    % ylim(rad2deg(lim_el));
    % xlabel('$\theta$');
    % ylabel('$\phi$');
    
    % annotation
    box on;
    ax=gca;
    set(ax,'Layer','Top');
    ax.FontSize=fontsize;
    ax.LineWidth=ax_lwidth;
    
    cbar=colorbar('SouthOutside');
    cbar.TickLabelInterpreter='latex';
    cbar.Label.Interpreter='latex';
    cbar.Label.String='$\Delta \mathrm{B}$ [mG]';
    cbar.Label.FontSize=fontsize;
    
    % save fig
    if do_save_figs
        savefigname=sprintf('fig_%s_%s',h.Name,getdatetimestr);
        fpath=fullfile(dir_save,savefigname);
        
        saveas(h,strcat(fpath,'.fig'),'fig');
        print(h,strcat(fpath,'.svg'),'-dsvg');
    end
end

%% Statistical uncertainty
if do_fit
    h=figure('Name','deltaB_unc','Units',f_units,'Position',f_pos,'Renderer',f_ren);
    
    [vazf,velf,deltaB_sef]=autofill_cent_symm(vaz,vel,deltaB_se_filt);
    tp=plotFlatMap(rad2deg(velf),rad2deg(vazf),1e3*deltaB_sef,'eckert4','texturemap');
    
    mu_se_dB=1e3*meanall(deltaB_sef,'omitnan');
    sig_se_dB=1e3*stdall(deltaB_sef,'omitnan');
    
    % annotation
    box on;
    ax=gca;
    set(ax,'Layer','Top');
    ax.FontSize=fontsize;
    ax.LineWidth=ax_lwidth;
    
    cbar=colorbar('SouthOutside');
    cbar.TickLabelInterpreter='latex';
    cbar.Label.Interpreter='latex';
    cbar.Label.String='$\mathrm{SE}(\Delta \mathrm{B})$ [mG]';
    cbar.Label.FontSize=fontsize;
    
    titlestr=sprintf('%s%0.1g(%0.1g)',...
        '$\overline{\mathrm{SE}}=$',mu_se_dB,sig_se_dB);  % label expparam
    title(titlestr);
    
    % save fig
    if do_save_figs
        savefigname=sprintf('fig_%s_%s',h.Name,getdatetimestr);
        fpath=fullfile(dir_save,savefigname);
        
        saveas(h,strcat(fpath,'.fig'),'fig');
        print(h,strcat(fpath,'.svg'),'-dsvg');
    end
end

%% VIS: tomography: equatorial
figname='dB_equatorial_tomography';
h=figure('Name',figname,'Units',f_units,'Position',f_pos,'Renderer',f_ren);

dB_eq=deltaB_filt(:,iel_0);      % delta magnetic field around equator
dBerr_eq=deltaB_se_filt(:,iel_0);      

hold on;
tp=ploterr(rad2deg(Vaz),1e3*dB_eq,[],1e3*dBerr_eq,'o','hhxy',0);
set(tp(1),'Marker',mark_typ{2},...
    'MarkerFaceColor',cl(2,:),'MarkerEdgeColor',c(2,:),...
    'DisplayName','');
set(tp(2),'Color',c(2,:));
pleg=tp(1);


% annotation
box on;
ax=gca;
set(ax,'Layer','Top');
ax.FontSize=fontsize;
ax.LineWidth=ax_lwidth;

title('Equatorial tomography $\phi=0$');

xlabel('Azimuthal angle $\theta$ [$^\circ$]');
ylabel('$\Delta B$ [mG]');

% axis tight;
xlim([0,180]);
ax.XTick=0:45:180;

% lgd=legend(pleg);

% save fig
if do_save_figs
    savefigname=sprintf('fig_%s_%s',figname,getdatetimestr);
    fpath=fullfile(dir_save,savefigname);
    
    saveas(h,strcat(fpath,'.fig'),'fig');
    print(h,strcat(fpath,'.svg'),'-dsvg');
end


%% integrated G
%% vis: G (integrated sum)
for ii=1:n_tau
    figname=sprintf('Gsum_%0.2fms',tau(ii));
    h=figure('Name',figname,'Units',f_units,'Position',f_pos,'Renderer',f_ren);
    
    [vazf,velf,tGf]=autofill_cent_symm(vaz,vel,squeeze(G_sum(ii,:,:)));
    tp=plotFlatMap(rad2deg(velf),rad2deg(vazf),tGf,'eckert4','texturemap');
    
    % annotation
    ax=gca;
    set(ax,'Layer','Top');
    ax.FontSize=fontsize;
    ax.LineWidth=ax_lwidth;
        
    titlestr=sprintf('%s%0.2f %s','$\tau=$',tau(ii),'ms');  % label expparam
    title(titlestr);
    
    cbar=colorbar('SouthOutside');
    cbar.TickLabelInterpreter='latex';
    cbar.Label.Interpreter='latex';
    cbar.Label.String='$\bar{g}^{(2)}$';
    cbar.Label.FontSize=fontsize;
    
    % save fig
    if do_save_figs
        savefigname=sprintf('fig_%s_%s',figname,getdatetimestr);
        fpath=fullfile(dir_save,savefigname);
        
        saveas(h,strcat(fpath,'.fig'),'fig');
        print(h,strcat(fpath,'.svg'),'-dsvg');
    end
end

%% vis: uncertainty of G 
for ii=1:n_tau    
    figname=sprintf('Gsum_unc_%0.2fms',tau(ii));
    h=figure('Name',figname,'Units',f_units,'Position',f_pos,'Renderer',f_ren);
    
    [vazf,velf,tGerrf]=autofill_cent_symm(vaz,vel,squeeze(G_sum_se(ii,:,:)));
    tp=plotFlatMap(rad2deg(velf),rad2deg(vazf),tGerrf,'eckert4','texturemap');
    
    mu_seG=meanall(G_sum_se(ii,:,:));
    sig_seG=stdall(G_sum_se(ii,:,:));
    
    % annotation
    ax=gca;
    set(ax,'Layer','Top');
    ax.FontSize=fontsize;
    ax.LineWidth=ax_lwidth;

    titlestr=sprintf('%s%0.2f %s; %s%0.1g(%0.1g)','$\tau=$',tau(ii),'ms',...
        '$\overline{\mathrm{SE}}=$',mu_seG,sig_seG);  % label expparam
    title(titlestr);

    cbar=colorbar('SouthOutside');
    cbar.TickLabelInterpreter='latex';
    cbar.Label.Interpreter='latex';
    cbar.Label.String='$\mathrm{SE}(\bar{g}^{(2)})$';
    cbar.Label.FontSize=fontsize;
    
    % save fig
    if do_save_figs
        savefigname=sprintf('fig_%s_%s',figname,getdatetimestr);
        fpath=fullfile(dir_save,savefigname);
        
        saveas(h,strcat(fpath,'.fig'),'fig');
        print(h,strcat(fpath,'.svg'),'-dsvg');
    end
end

%% statistics
% TODO:
%   azim-elev is latlon grid which gives a *non-uniform* sampling
%

for ii=1:n_tau
    figname=sprintf('hist_Gsum_%0.2fms',tau(ii));
    h=figure('Name',figname,'Units',f_units,'Position',f_pos,'Renderer',f_ren);
    
    thist=histogram(G_sum(ii,:,:,:),'Normalization','pdf');
%     thist.LineWidth=ax_lwidth;
    
    hold on;
    % Gaussian fit
    tmu=meanall(G_sum(ii,:,:,:));
    tsigma=stdall(G_sum(ii,:,:,:));
    ax_xlim=xlim;
    x=linspace(ax_xlim(1),ax_xlim(2));
    y=exp(-(x-tmu).^2./(2*tsigma^2))./(tsigma*sqrt(2*pi));
    str_gauss=sprintf('%s(%0.2g,%0.1g)','$(\mu,\sigma)=$',tmu,tsigma);
    p_gauss=plot(x,y,'LineWidth',line_wid,...
        'DisplayName',str_gauss);
    
    % annotation
    ax=gca;
    set(ax,'Layer','Top');
    ax.FontSize=fontsize;
%     ax.LineWidth=ax_lwidth;
    
    titlestr=sprintf('%s%0.2f %s','$\tau=$',tau(ii),'ms');  % label expparam
    title(titlestr);
    
    xlabel('$\bar{g}^{(2)}$');
    ylabel('Norm frac lat-lon zones');
    
    lgd=legend(p_gauss);
    
    % save fig
    if do_save_figs
        savefigname=sprintf('fig_%s_%s',figname,getdatetimestr);
        fpath=fullfile(dir_save,savefigname);
        
        saveas(h,strcat(fpath,'.fig'),'fig');
        print(h,strcat(fpath,'.svg'),'-dsvg');
    end
end

%% End of script
toc