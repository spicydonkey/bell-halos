%% Time evolution of Psi+ halo
% DKS
% 2018-10-25

tic

%% CONFIGS
% data file
fdata='C:\Users\HE BEC\Dropbox\PhD\projects\halo_metrology\analysis\exp4_tevo\exp4_20181029.mat';


% halo & collision  ---------------------------------------
v_sep=120e-3;       % pair separation velocity [m/s]

% Ramsey (3ms)
d_sep = 3e-3 * v_sep;      % 3ms expansion
sig_psf_r = 35e-6;      % rms width of point-spread function ~ 35um (BEC) (see supplementary in paper)

sig_psf_beta=sig_psf_r/(d_sep/2);   % PSF rms width in angle (rad)

% gradiometry (~1.5 ms)
sig_psf_gradiometry=sig_psf_beta*2;     % ~twice uncertainty


% grid/bins ------------------------------------------------------
% alpha=pi/8;        % cone half-angle
alpha=sig_psf_gradiometry;
lim_az=[0,pi];              % limit of azim angle (exc +pi)
phi_max=pi/4;
lim_el=[-phi_max,phi_max];

n_az=40;                	% equispaced bins
n_el=20;

az_disp=deg2rad([0,45,90]);     % azim sections (great circles) to display


% g2 -------------------------------------------------------------
n_dk=7;
lim_dk=[-0.2,0.2];

% bootstrapping
bs_frac=0.1;      
bs_nrep=20;


% Physical constants ---------------------------------------------
Cphys=physConsts;
C_gymag=Cphys.He_gymag;         % He* gyromagnetic ratio (gamma) [Hz/G]


% vis ----------------------------------------------------------------
config_fig = loadFigureConfig;

f_units='normalized';
f_pos=[0.2,0.2,0.2,0.3];
f_pos_wide=[0.2,0.2,0.25,0.3];
f_ren='painters';

[c,cl,cd]=palette(3);
line_sty={'-','--',':','-.'};
mark_typ={'+','o','^','d'};
mark_siz=7;
line_wid=1.5;
fontsize=12;
ax_lwidth=1.2;

% theme color
ctheme=parula(5);
ctheme=ctheme(1:end-1,:);
cltheme=colshades(ctheme);

idx_col=2;     % mid-color 

% marker colors
c_loc=zeros(length(iaz_disp),3);        % *BLACK* markers (simple)
cl_loc=[1,1,0]'*[1,1,1];

%% load data
load(fdata);

%% transform to RH coords
% SAVE ORIGINAL 
if ~exist('k_tau_orig','var')
    k_tau_orig=k_tau;
end

k_tau=cellfun(@(C) cellfun(@(x) tzxy2RHtzxy2(x),C,'UniformOutput',false),...
    k_tau_orig,'UniformOutput',false);      % EXP-coord sys (z against g)


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

%% vis: Parity - ALL
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
xlabel('$\tau~(\textrm{ms})$');
ylabel('Parity $\bar{\mathcal{B}}_{\pi/2}$');
ax.FontSize=fontsize;
ax.LineWidth=ax_lwidth;
ylim(1.5*[-1,1]);


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
        
        %%% fitted model
%         % 1
%         tpf=plot(tau0_fit+T_so,B0_fit{iaz,jj},'-',...
%             'LineWidth',line_wid,'Color',cc(jj,:));
%         uistack(tpf,'bottom');
        
        % 2
        tpf=plot(t_fit,B0_fit2{iaz,jj},'-',...
            'LineWidth',line_wid,'Color',cc(jj,:));
        uistack(tpf,'bottom');
        
    end
    
    % annotation
    box on;
    ax=gca;
    set(ax,'Layer','Top');
    xlabel('$\tau~(\textrm{ms})$');
    ylabel('Parity $\bar{\mathcal{B}}_{\pi/2}$');
    ax.FontSize=fontsize;
    ax.LineWidth=1.2;
    ylim([-1.2,1.2]);
    titlestr=sprintf('%s %0.0f','$\theta=$',rad2deg(taz));
    title(titlestr);
    lgd=legend(pleg,'Location','EastOutside');
    title(lgd,'Latitude $\phi$ (deg)');
end


%% vis: Parity distribution (3D)
H=[];
for ii=1:n_tau
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
    ax.CLim=1.5*[-1,1];
    cbar.Limits=ax.CLim;
    cbar.TickLabelInterpreter='latex';
    cbar.Label.Interpreter='latex';
    cbar.Label.String='Parity $\bar{\mathcal{B}}_{\pi/2}$';
    cbar.Label.FontSize=fontsize;
end

%% vis: Parity distribution (2D)
H=[];
for ii=1:n_tau
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
    ax.CLim=1.5*[-1,1];
    cbar.Limits=ax.CLim;
    cbar.TickLabelInterpreter='latex';
    cbar.Label.Interpreter='latex';
    cbar.Label.String='Parity $\bar{\mathcal{B}}_{\pi/2}$';
    cbar.Label.FontSize=fontsize;
end


%% Gradiometry Analysis ================================================

%% Model 1: time-independent asymmetry of Bell triplet
% % B_pi/2 = cos(hbar*gamma* deltaB * t)
% 
% % set up model and solver
% mdl_tevo.mdl='y~cos(om*x)';
% mdl_tevo.cname={'om'};
% mdl_tevo.fopt=statset('TolFun',10^-10,'TolX',10^-10,'MaxIter',10^6,'UseParallel',0);
% par0=1;
% 
% % preproc
% T_so=0.8;           % MODEL PARAM: asymmetry "switch-on" time = 0.8 ms
% tau0=tau-T_so;      % time evolution since switch-on
% 
% % fit
% mdl_tevo.fit=cell(n_az,n_el);
% for ii=1:n_zone
%     [iaz,iel]=ind2sub([n_az,n_el],ii);
%     
%     tB0=B0(:,iaz,iel);
%     mdl_tevo.fit{iaz,iel}=fitnlm(tau0,tB0,mdl_tevo.mdl,par0,...
%         'CoefficientNames',mdl_tevo.cname,'Options',mdl_tevo.fopt);
% end
% 
% % eval fitted model
% tau0_fit=linspace(0,5,1e3);     % tau vals to eval fit
% mdl_tevo.fit_par=cellfun(@(m) m.Coefficients.Estimate,mdl_tevo.fit);
% mdl_tevo.fit_par_se=cellfun(@(m) m.Coefficients.SE,mdl_tevo.fit);
% B0_fit=cellfun(@(f) feval(f,tau0_fit),mdl_tevo.fit,'UniformOutput',false);
% 
% % diff B field 
% om_fit=1e3*mdl_tevo.fit_par;        % fitted omega [rad/s]
% om_se_fit=1e3*mdl_tevo.fit_par_se;
% deltaB=om_fit/(2*pi*C_gymag);        % diff in B-field strength [G]
% deltaB_se=om_se_fit/(2*pi*C_gymag);   % standard error (fit)


%% Model 2: constant B-field gradient + variable creation time
% MODEL:
%      B(x) = B0 + dBdx * x    --- const grad mag-field
%      x(t) = v * (t - t0)     --- const velocity; created at t0
% DeltaB(t) = B(xA) - B(xB) = dBdx * (2*v*(t-t0))    --- diff B-field
%
%%%
% parity vs. phase
%   parity = cos( 2*PHI(t) )
%  where 
%   PHI(t) = C_gymag \int DeltaB(t) dt
%          = C_gymag * dBdx * v * (t-t0)^2
%          = beta * (t-t0)^2
%
% CHECK ABOVE

%-------------------------------------------------------------------
% t0 is constrained to [0,0.4] ms
% mdl_tevo2.mdl='y~cos(2*beta*(x-0.4*(tanh(phi)+1)/2)^2)';
% mdl_tevo2.cname={'beta','phi'};     
% mdl_tevo2.fopt=statset('TolFun',10^-10,'TolX',10^-10,'MaxIter',10^6,'UseParallel',0);
% mdl_tevo2.par0=[0.15,0];

%--------------------------------------------------------------
% t0 is fully constrained
mdl_tevo2.mdl='y~cos(2*beta*(x)^2)';       
mdl_tevo2.cname={'beta'};     
mdl_tevo2.fopt=statset('TolFun',10^-10,'TolX',10^-10,'MaxIter',10^6,'UseParallel',0);
mdl_tevo2.par0=[0.15];

% ---------------------------------------------------------------
% TRY: parity amplitude is free
% mdl_tevo2.mdl='y~amp*cos(2*beta*(x-0.4)^2)';       % constrain t0
% mdl_tevo2.cname={'beta','amp'};     
% mdl_tevo2.fopt=statset('TolFun',10^-10,'TolX',10^-10,'MaxIter',10^6,'UseParallel',0);
% mdl_tevo2.par0=[0.15,1];

% mdl_tevo2.mdl='y~amp*cos(2*beta*(x-0.8*(tanh(phi)+1)/2)^2)';       % constrain t0 to [0,0.8]
% mdl_tevo2.cname={'beta','phi','amp'};     
% mdl_tevo2.fopt=statset('TolFun',10^-10,'TolX',10^-10,'MaxIter',10^6,'UseParallel',0);
% mdl_tevo2.par0=[0.15,0,1];


%% FIT: Model 2
mdl_tevo2.fit=cell(n_az,n_el);
for ii=1:n_zone
    [iaz,iel]=ind2sub([n_az,n_el],ii);
    
    tB0=B0(:,iaz,iel);
    mdl_tevo2.fit{iaz,iel}=fitnlm(tau,tB0,mdl_tevo2.mdl,mdl_tevo2.par0,...
        'CoefficientNames',mdl_tevo2.cname,'Options',mdl_tevo2.fopt);
end

% eval fitted model
t_fit=linspace(0,5,1e3);     % time since collision pulse to eval fit
mdl_tevo2.fit_par=arrayfun(@(I) cellfun(@(m) m.Coefficients.Estimate(I),mdl_tevo2.fit),1:numel(mdl_tevo2.cname),'UniformOutput',false);
mdl_tevo2.fit_par_se=arrayfun(@(I) cellfun(@(m) m.Coefficients.SE(I),mdl_tevo2.fit),1:numel(mdl_tevo2.cname),'UniformOutput',false);
B0_fit2=cellfun(@(f) feval(f,t_fit),mdl_tevo2.fit,'UniformOutput',false);


%%% Magnetic field gradient: dBdx -----------------------------------------
beta=abs(mdl_tevo2.fit_par{1});             % get fit params and ensure sign is positive
beta_se=abs(mdl_tevo2.fit_par_se{1});

% far-field
dBdx=1e6*beta/(v_sep*C_gymag/2);            % [G/m] (factor to convert time-scaling ms^-2 --> s^-2)
dBdx_se=1e6*beta_se/(v_sep*C_gymag/2);      

% fill by inversion symmetry
[gaz,gel,dBdx]=autofill_cent_symm(vaz,vel,dBdx);
[~,~,dBdx_se]=autofill_cent_symm(vaz,vel,dBdx_se);

% define 2pi-wrap filled az-el vectors
az=gaz(:,1)';
iaz_0=find(az==0);
el=gel(1,:);

% % outlier to NaN (ones that have 1e4 mean errors)
% b_outlier=isoutlier(dBdx_se_ff);
% dBdx_ff(b_outlier)=NaN;
% dBdx_se_ff(b_outlier)=NaN;

% % interrogation time
% % METHOD: estimated at interrogation time
% %   Gaussian convolution: see SMs for determination of length scale
% %       ERROR by conv of variance
% dBdx_r=gaussfilt_sph(dBdx_ff,gaz,gel,gaz,gel,sig_psf_beta);      
% dBdxerr_r=sqrt(gaussfilt_sph(dBdx_se_ff.^2,gaz,gel,gaz,gel,sig_psf_beta));


%% VIS: fit params diagnostic
h=figure('Name','fit_param_diagnostic');
n_fit_param=numel(mdl_tevo2.fit_par);
for ii=1:n_fit_param
    % estimate
    ax=subplot(2,n_fit_param,sub2ind([n_fit_param,2],ii,1));
    imagesc(rad2deg(az),rad2deg(el),mdl_tevo2.fit_par{ii}');
    cbar=colorbar();
    title(cbar,mdl_tevo2.cname{ii});
    ax.YDir='normal';
    title(mdl_tevo2.cname{ii});
    
    % error
    ax=subplot(2,n_fit_param,sub2ind([n_fit_param,2],ii,2));
    imagesc(rad2deg(az),rad2deg(el),mdl_tevo2.fit_par_se{ii}');
    cbar=colorbar();
    title(cbar,'SE');
    ax.YDir='normal';
    
    %     xlabel('$\theta$');
    %     ylabel('$\phi$');
end


%% STAT: regions to display parity signal
% max, min and middle
loc_disp_2pi=[];    % row-array of (Iaz,Iel)-location to display (in 2*pi wrapped)
azel_disp=[];   % (az,el)

%%% theta,phi=(0,0) ----------------------------------------------
% loc_disp_2pi(end+1,:)=[iaz_0,iel_0];
% azel_disp(end+1,:)=[az(iaz_0),el(iel_0)];

%%% max -----------------------------------------------------------
dBdx_max=max(dBdx(:));
[iaz_max,iel_max]=find(dBdx==dBdx_max);
iaz_max=iaz_max(end);     % due to symmetry there can be multiple
iel_max=iel_max(end);
loc_disp_2pi(end+1,:)=[iaz_max,iel_max];
azel_disp(end+1,:)=[az(iaz_max),el(iel_max)];

%%% min --------------------------------------------------------------
dBdx_min=min(dBdx(:));
[iaz_min,iel_min]=find(dBdx==dBdx_min);
iaz_min=iaz_min(1);
iel_min=iel_min(1);
loc_disp_2pi(end+1,:)=[iaz_min,iel_min];
azel_disp(end+1,:)=[az(iaz_min),el(iel_min)];

%%% middle (in-between min-max) ----------------------------------------
iaz_mid=round(0.5*(iaz_max+iaz_min));
iel_mid=round(0.5*(iel_max+iel_min));
loc_disp_2pi(end+1,:)=[iaz_mid,iel_mid];
azel_disp(end+1,:)=[az(iaz_mid),el(iel_mid)];

% summary
n_loc_disp=size(loc_disp_2pi,1);        % num of locations to display

%%% Find location in the original hemisphere Vaz,Vel -------------------
loc_disp_orig=NaN(size(loc_disp_2pi));
azel_disp_orig=azel_disp;
for ii=1:n_loc_disp
    taz=azel_disp(ii,1);
    tel=azel_disp(ii,2);
    
    % check if in orig hemisphere
    if ~all(any(Vaz==taz) && any(Vel==tel))
        % spherical inversion
        azel_disp_orig(ii,:)=[wrapTo2Pi(taz+pi),-tel];
    end
        
    % get index location
    loc_disp_orig(ii,1)=find(Vaz==azel_disp_orig(ii,1));
    loc_disp_orig(ii,2)=find(Vel==azel_disp_orig(ii,2));
end


%% vis: Parity evolution fit (publication)
figname='B0_tevo_eqt';
% figure('Name',figname,'Units',f_units,'Position',[0.2,0.2,0.5,0.2],'Renderer',f_ren);
h=figure('Name',figname,'Units','centimeters','Position',[0,0,8.6,3.2],'Renderer',f_ren);
hold on;

% disp_iaz=iaz_disp;

pleg=NaN(n_loc_disp,1);
for ii=1:n_loc_disp
    iazel=loc_disp_orig(ii,:);
    tazel=azel_disp_orig(ii,:);
    
%     tiaz=disp_iaz(ii);
% tp=ploterr(tau,squeeze(B0(:,tiaz,iel_0)),[],squeeze(B0_bs_se(:,tiaz,iel_0)),'o','hhxy',0);

    tp=ploterr(tau,squeeze(B0(:,iazel(1),iazel(2))),[],squeeze(B0_bs_se(:,iazel(1),iazel(2))),'o','hhxy',0);
    set(tp(1),'Marker',mark_typ{ii},'MarkerSize',4.5,...
        'MarkerFaceColor',cl_loc(ii,:),'Color',c_loc(ii,:),'DisplayName',num2str(rad2deg(Vaz(tiaz)),2));
    set(tp(2),'Color',c_loc(ii,:));
    pleg(ii)=tp(1);
    
    %%% fitted model
%     %%%% MODEL1
%     % stationary before tau0 (Psi+ stationary --> parity=1)
%     tpf=plot([0,T_so],[1,1],...
%         line_sty{ii},'LineWidth',config_fig.line_wid,'Color',c_loc(ii,:));
%     uistack(tpf,'bottom');
%     
%     % dynamics after turn-on
%     tpf=plot(tau0_fit+T_so,B0_fit{tiaz,iel_0},...
%         line_sty{ii},'LineWidth',config_fig.line_wid,'Color',c_loc(ii,:));
%     uistack(tpf,'bottom');

    %%%% MODEL2
    tpf=plot(t_fit,B0_fit2{iazel(1),iazel(2)},...
    line_sty{ii},'LineWidth',config_fig.line_wid,'Color',c_loc(ii,:));
    uistack(tpf,'bottom');
end

% annotation
box on;
ax=gca;
set(ax,'Layer','Top');
xlabel('$\tau~(\textrm{ms})$');
% ylabel('Parity $\bar{\mathcal{B}}_{\pi/2}$');
ylabel('parity');

ax.FontSize=config_fig.ax_fontsize;
ax.LineWidth=config_fig.ax_lwid;

xlim([0.7,1.8]);
ylim(1.5*[-1,1]);
% axis auto;

% lgd=legend(pleg,'Location','SouthWest');
% title(lgd,'Azimuth $\theta$ (deg)');


%% vis: Model 1: deltaB (asymmetry measure) around halo: t-indep model (3D SPH)
% h=figure('Name','deltaB_sphdist_3d','Units',f_units,'Position',f_pos,'Renderer',f_ren);
% 
% % tp=plot_sph_surf(vaz,vel,1e3*deltaB);
% [vazf,velf,deltaBf]=autofill_cent_symm(vaz,vel,deltaB);
% tp=plot_sph_surf(vazf,velf,1e3*deltaBf);
% 
% % annotation
% ax=gca;
% ax.FontSize=fontsize;
% 
% axis on;
% box on;
% xlabel('$k_x$');
% ylabel('$k_y$');
% zlabel('$k_z$');
% xlim([-1,1]);
% ylim([-1,1]);
% zlim([-1,1]);
% 
% cbar=colorbar();
% cbar.TickLabelInterpreter='latex';
% cbar.Label.Interpreter='latex';
% cbar.Label.String='$\Delta \mathrm{B}$ (mG)';
% cbar.Label.FontSize=fontsize;

%% vis: dBdx around halo: Model2 (3D SPH)
% h=figure('Name','dBdx_sphdist_3d','Units',f_units,'Position',f_pos,'Renderer',f_ren);
% 
% % [vazf,velf,dBdxf]=autofill_cent_symm(vaz,vel,dBdx_ff);
% % tp=plot_sph_surf(vazf,velf,dBdxf);
% tp=plot_sph_surf(gaz,gel,dBdx);
% 
% 
% % annotation
% ax=gca;
% ax.FontSize=fontsize;
% 
% axis on;
% box on;
% xlabel('$k_x$');
% ylabel('$k_y$');
% zlabel('$k_z$');
% xlim([-1,1]);
% ylim([-1,1]);
% zlim([-1,1]);
% 
% cbar=colorbar();
% cbar.TickLabelInterpreter='latex';
% cbar.Label.Interpreter='latex';
% cbar.Label.String='$d\mathrm{B}/dx$ (G/m)';
% cbar.Label.FontSize=fontsize;

%% vis: deltaB (asymmetry measure) around halo: t-indep model (2D PROJ MAP)
% h=figure('Name','deltaB_sphdist_2d','Units',f_units,'Position',f_pos,'Renderer',f_ren);
% 
% % tp=plotFlatMap(rad2deg(vel),rad2deg(vaz),1e3*deltaB,'rect','texturemap');
% tp=plotFlatMap(rad2deg(velf),rad2deg(vazf),1e3*deltaBf,'eckert4','texturemap');
% 
% % all around the halo (redundant since inversion symmetry)
% % tp=plotFlatMapWrappedRad(vazf,velf,1e3*deltaBf,'rect','texturemap');
% 
% % annotation
% ax=gca;
% ax.FontSize=fontsize;
% 
% cbar=colorbar('SouthOutside');
% cbar.TickLabelInterpreter='latex';
% cbar.Label.Interpreter='latex';
% cbar.Label.String='$\Delta \mathrm{B}$ (mG)';
% cbar.Label.FontSize=fontsize;


%% vis: dBdx - rectangular (publication) (@ff)
figname='dBdx_rectmap_ff';
h=figure('Name',figname,'Units','centimeters','Position',[0,0,8.6,4.5],'Renderer',f_ren);

tp=plotFlatMapWrappedRad(gaz,gel,dBdx,'rect','texturemap');
% tp=plotFlatMapWrappedRad(gaz,gel,dBdx_se_ff,'rect','texturemap');

% % AZIM in [0,pi]
% vaz_pi=vaz;
% vaz_pi(end+1,:)=vaz_pi(1,:)+pi;
% vel_pi=vel;
% vel_pi(end+1,:)=vel_pi(1,:);
% dBdx_pi=dBdx;
% dBdx_pi(end+1,:)=fliplr(dBdx_pi(1,:));
% tp=plotFlatMap(rad2deg(vel_pi),rad2deg(vaz_pi),dBdx_pi,'rect','texturemap');
% tp=plotFlatMap(rad2deg(vel),rad2deg(vaz),dBdx_r,'rect','texturemap');

% label ROI
hold on;
for ii=1:n_loc_disp
    iazel=loc_disp_2pi(ii,:);
    tazel=azel_disp(ii,:);
    
%     tiaz=disp_iaz(ii);
    tp=plot(rad2deg(tazel(1)),rad2deg(tazel(2)),'Marker',mark_typ{ii},...
        'MarkerEdgeColor',c_loc(ii,:),'MarkerFaceColor',cl_loc(ii,:),...
        'MarkerSize',config_fig.mark_siz);
end

% annotation
ax=gca;
set(ax,'Layer','Top');
box on;
% grid on;
ax.FontSize=config_fig.ax_fontsize;
ax.LineWidth=config_fig.ax_lwid;

% axis tight;
xlim([-180,180]);
ylim([-90,90]);
xticks(-180:90:180);
yticks(-90:45:90);

xlabel('$\theta$ (deg)');
ylabel('$\phi$ (deg)');

colormap('parula');
cbar=colorbar('eastoutside');
clim_original=cbar.Limits;
cbar.Limits=[0,clim_original(2)];
cbar.TickLabelInterpreter='latex';
cbar.Label.Interpreter='latex';
cbar.Label.String='$d\mathrm{B}/dx$ (G/m)';
cbar.Label.FontSize=config_fig.ax_fontsize;
cbar.FontSize=config_fig.ax_fontsize;

% change colorbar width
pos_ax=get(gca,'Position');
pos_cbar=get(cbar,'Position');
pos_cbar(3)=0.03;
set(cbar,'Position',pos_cbar);
set(gca,'Position',pos_ax);


% hatch-out truncated region --------------------------
% % hack: two separate hatched regions
% hpatch_trunc(1)=patch('XData',180*[-1,1,1,-1],...
%     'YData',[-90,-90,-45,-45],...
%     'FaceColor','none','EdgeColor','none');
% H_trunc(1)=hatchfill2(hpatch_trunc(1),'single','HatchAngle',45,'HatchDensity',20,...
%     'HatchColor','k','HatchLineWidth',config_fig.line_wid);
% 
% hpatch_trunc(2)=patch('XData',180*[-1,1,1,-1],...
%     'YData',[45,45,90,90],...
%     'FaceColor','none','EdgeColor','none');
% H_trunc(2)=hatchfill2(hpatch_trunc(2),'single','HatchAngle',45,'HatchDensity',20,...
%     'HatchColor','k','HatchLineWidth',config_fig.line_wid);

% TODO - vecrast (the whole axes as bottom layer)
% hpatch_trunc=patch('XData',[ax.XLim(1),ax.XLim(2),ax.XLim(2),ax.XLim(1)],...
%     'YData',[ax.YLim(1),ax.YLim(1),ax.YLim(2),ax.YLim(2)],...
%     'FaceColor','none','EdgeColor','none');
% uistack(hpatch_trunc,'bottom');
% H_trunc=hatchfill2(hpatch_trunc,'single','HatchAngle',45,'HatchDensity',20,...
%     'HatchColor','k','HatchLineWidth',config_fig.line_wid);


%---------------------------------------------
% NOTE: save with vecrast
% e.g.
% vecrast(h,strcat('dB_dist_',getdatetimestr),600,'bottom','pdf')
%---------------------------------------------


%% vis: Model 1: deltaB - rectangular (publication)
% figname='deltaB_rectmap';
% h=figure('Name',figname,'Units','centimeters','Position',[0,0,8.6,4.5],'Renderer',f_ren);
% 
% % AZIM in [0,pi]
% vaz_pi=vaz;
% vaz_pi(end+1,:)=vaz_pi(1,:)+pi;
% vel_pi=vel;
% vel_pi(end+1,:)=vel_pi(1,:);
% deltaB_pi=deltaB;
% deltaB_pi(end+1,:)=fliplr(deltaB_pi(1,:));
% tp=plotFlatMap(rad2deg(vel_pi),rad2deg(vaz_pi),1e3*deltaB_pi,'rect','texturemap');
% % tp=plotFlatMap(rad2deg(vel_pi),rad2deg(vaz_pi),1e3*deltaB_pi,'eckert4','texturemap');
% 
% % label ROI
% hold on;
% for ii=1:numel(disp_iaz)
%     tiaz=disp_iaz(ii);
%     tp=plot(rad2deg(Vaz(tiaz)),rad2deg(Vel(iel_0)),'Marker',mark_typ{ii},...
%         'MarkerEdgeColor',c_loc(ii,:),'MarkerFaceColor',cl_loc(ii,:),...
%         'MarkerSize',config_fig.mark_siz);
% end
% 
% % annotation
% ax=gca;
% set(ax,'Layer','Top');
% box on;
% % grid on;
% ax.FontSize=config_fig.ax_fontsize;
% ax.LineWidth=config_fig.ax_lwid;
% 
% axis tight;
% % xlim([-180,180]);
% xticks(-180:90:180);
% yticks(-90:45:90);
% 
% xlabel('$\theta$ (deg)');
% ylabel('$\phi$ (deg)');
% 
% % ax.DataAspectRatio=[1 0.5 1];       % adjust aspect ratio
% 
% cbar=colorbar('eastoutside');
% cbar.TickLabelInterpreter='latex';
% cbar.Label.Interpreter='latex';
% cbar.Label.String='$\Delta \mathrm{B}$ (mG)';
% cbar.Label.FontSize=config_fig.ax_fontsize;
% cbar.FontSize=config_fig.ax_fontsize;
% 
% % change colorbar width
% pos_ax=get(gca,'Position');
% pos_cbar=get(cbar,'Position');
% pos_cbar(3)=0.03;
% set(cbar,'Position',pos_cbar);
% set(gca,'Position',pos_ax);
% 
% %---------------------------------------------
% % NOTE: save with vecrast
% % e.g.
% % vecrast(h,strcat('dB_dist_',getdatetimestr),600,'bottom','pdf')
% %---------------------------------------------


%% VIS: Model 1: Diff B: Equatorial tomography (publication)
% %%% calibration
% % nuller characterisation at ±y measurement
% dBdy=0.4;           % mG/mm
% dBdy_ferr=0.6;      % fractional error
% th_exp=pi/2;        % deltaY - measurement direction
% % th_err_exp        % TODO: uncertainty of orientation
% 
% 
% %%% configs
% p_xlim=[0,180];     % periodic BC
% 
% % ctheme=magma(5);
% ctheme=parula(5);
% ctheme=ctheme(1:end-1,:);
% cltheme=colshades(ctheme);
% 
% idx_col=2;     % mid-color 
% 
% %%% DATA
% dB_eq=deltaB(:,iel_0);      % delta magnetic field around equator
% dBerr_eq=deltaB_se(:,iel_0);      
% 
% % fit
% sine_mdl='y~c+amp*cos(2*x+phi)';      % periodic boundary condition fixes OMEGA
% sine_cname={'amp','phi','c'};
% sine_par0=[0.05,0,0.1];
% 
% fit_dB_eq=fitnlm(Vaz,1e3*dB_eq,sine_mdl,sine_par0,'CoefficientNames',sine_cname);
% xx=linspace(0,pi,1e3);
% yy=feval(fit_dB_eq,xx);
% 
% %%% PLOT
% figname='dB_equatorial_tomography';
% % h=figure('Name',figname,'Units',f_units,'Position',[0.2,0.2,0.5,0.2],'Renderer',f_ren);
% h=figure('Name',figname,'Units','centimeters','Position',[0,0,8.6,3.2],'Renderer',f_ren);
% 
% hold on;
% 
% % data
% %%% scatter with errbars
% % skip ROI data
% not_roi = true(size(Vaz));
% not_roi(iaz_disp)=false;
% 
% tp=ploterr(rad2deg(Vaz(not_roi)),1e3*dB_eq(not_roi),[],1e3*dBerr_eq(not_roi),'.','hhxy',0);
% set(tp(1),'Marker','*','MarkerSize',config_fig.mark_siz,...
%     'MarkerFaceColor',[1,1,1],'MarkerEdgeColor',ctheme(idx_col,:));
% set(tp(2),'Color',ctheme(idx_col,:));
% pleg=tp(1);
% 
% % %%% SHADED ERR BAR
% % tp=shadedErrorBar(rad2deg(Vaz),1e3*dB_eq,1e3*dBerr_eq,'r');
% % tp.mainLine.Color='none';  %ctheme(idx_col,:);
% % tp.mainLine.LineWidth=config_fig.ax_lwid;
% % tp.patch.FaceColor=ctheme(idx_col,:);       %cltheme(idx_col,:);
% % tp.patch.FaceAlpha=0.33;
% % tp.edge(1).Visible='off';
% % tp.edge(2).Visible='off';
% 
% % fit
% pfit=plot(rad2deg(xx),yy,'LineStyle','-','Color',ctheme(idx_col,:),...
%     'LineWidth',config_fig.ax_lwid);
% uistack(pfit,'bottom');
% 
% % independent measurement (calibration)
% d_sep_calib=0.15;             % halo effective Dia. (separation) for diff B (mm)
% dB_exp=d_sep_calib*dBdy;      % diff B (mG)
% 
% tp_exp=ploterr(rad2deg(th_exp),dB_exp,[],dB_exp*dBdy_ferr,'.','hhxy',0);
% set(tp_exp(1),'Marker','none');
% set(tp_exp(2),'Color',0.8*ones(1,3),'LineWidth',5);
% arrayfun(@(p) uistack(p,'bottom'),tp_exp);
% 
% %%% ROI
% for ii=1:numel(disp_iaz)
%     clearvars tp;
%     
%     iazel=[disp_iaz(ii),iel_0];
%     taz_disp=Vaz(iazel(1));    
%     tp=ploterr(rad2deg(taz_disp),1e3*dB_eq(iazel(1)),[],1e3*dBerr_eq(iazel(1)),...
%         mark_typ{ii},'hhxy',1);
%     
%     set(tp(1),'MarkerEdgeColor',c_loc(ii,:),...
%         'MarkerFaceColor',cl_loc(ii,:),...
%         'MarkerSize',config_fig.mark_siz,...
%         'DisplayName',num2str(rad2deg([taz_disp,Vel(iel_0)])));
%     set(tp(2),'Color',c_loc(ii,:));
% end
% 
% 
% %%% annotation
% box on;
% ax=gca;
% set(ax,'Layer','Top');
% ax.FontSize=config_fig.ax_fontsize;
% ax.LineWidth=config_fig.ax_lwid;
% 
% xlabel('$\theta$ (deg)');
% ylabel('$\Delta \mathrm{B}$ (mG)');
% 
% ylim_0=ax.YLim;
% ylim([0,ylim_0(2)]);
% xlim(p_xlim);
% ax.XTick=0:45:180;
% 
% % lgd=legend(pleg);


%% VIS: dBdx Equatorial tomography (publication)
%%% calibration
% % nuller characterisation at ±y measurement
% dBdy=0.4;           % mG/mm
% dBdy_ferr=0.6;      % fractional error
% th_exp=pi/2;     % deltaY - measurement direction
% % th_err_exp    % TODO: uncertainty of orientation


%%% configs
p_xlim=[0,180];     % periodic BC
idx_col=2;     % mid-color 

%%% DATA
dBdx_eq=dBdx(:,iel_0);      % delta magnetic field around equator
dBdxerr_eq=dBdx_se(:,iel_0);      

% % fit
% sine_mdl='y~c+amp*cos(2*x+phi)';      % periodic boundary condition fixes OMEGA
% sine_cname={'amp','phi','c'};
% sine_par0=[0.05,0,0.1];
% 
% fit_dB_eq=fitnlm(Vaz,1e3*dB_eq,sine_mdl,sine_par0,'CoefficientNames',sine_cname);
% xx=linspace(0,pi,1e3);
% yy=feval(fit_dB_eq,xx);

%%% PLOT
figname='dBdx_equatorial_tomography';
h=figure('Name',figname,'Units','centimeters','Position',[0,0,8.6,3.2],'Renderer',f_ren);
hold on;

% data-------------------------------------------------------------------
% Ramsey estimate
ramsey_fname='C:\Users\HE BEC\Documents\MATLAB\bell-halos\analysis\ramsey_mag_tomo\out\out_20190220_142302.mat';

vars_to_load={'az','el','dBdx_eq','dBdx_eq_se'};
S_ramsey=load(ramsey_fname,vars_to_load{:});

% SHADED ERR BAR
tp=shadedErrorBar(rad2deg(S_ramsey.az),S_ramsey.dBdx_eq,S_ramsey.dBdx_eq_se,'k');
tp.mainLine.Color=0.3*ones(1,3);    %  'k';    % 'none';
tp.mainLine.LineWidth=1;
tp.mainLine.LineStyle='--';
% tp.patch.FaceColor=ctheme(idx_col,:);       %cltheme(idx_col,:);
% tp.patch.FaceAlpha=0.33;
tp.edge(1).Visible='off';
tp.edge(2).Visible='off';


% Ent-based Gradiometry
% % skip ROI data
% not_roi = true(size(az));
% not_roi(iaz_0+iaz_disp-1)=false;

% SHADED ERR BAR
tp=shadedErrorBar(rad2deg(az),dBdx_eq,dBdxerr_eq,'r');
tp.mainLine.Color=ctheme(idx_col,:);    % 'none';
tp.mainLine.LineWidth=1;
tp.patch.FaceColor=ctheme(idx_col,:);       %cltheme(idx_col,:);
tp.patch.FaceAlpha=0.33;
tp.edge(1).Visible='off';
tp.edge(2).Visible='off';

% %%% ROI
% for ii=1:numel(disp_iaz)
%     clearvars tp;
%     
%     iazel=[disp_iaz(ii),iel_0];
%     taz_disp=Vaz(iazel(1));    
%     tp=ploterr(rad2deg(taz_disp),dBdx_eq(iazel(1)),[],dBdxerr_eq(iazel(1)),...
%         mark_typ{ii},'hhxy',1);
%     
%     set(tp(1),'MarkerEdgeColor',c_loc(ii,:),...
%         'MarkerFaceColor',cl_loc(ii,:),...
%         'MarkerSize',config_fig.mark_siz,...
%         'DisplayName',num2str(rad2deg([taz_disp,Vel(iel_0)])));
%     set(tp(2),'Color',c_loc(ii,:));
% end

%%% annotation
box on;
ax=gca;
set(ax,'Layer','Top');
ax.FontSize=config_fig.ax_fontsize;
ax.LineWidth=config_fig.ax_lwid;

xlabel('$\theta$ (deg)');
ylabel('$d\mathrm{B}/dx$ (G/m)');

axis tight;
ylim_0=ax.YLim;
ylim([0,ylim_0(2)]);
xlim(p_xlim);
ax.XTick=0:45:180;

% lgd=legend(pleg);

%%% Zoomed - INSET ------------------------------------------------
ax2=axes('Position',[.7 .5 .2 .4]);      % inset axes 'YAxisLocation','right'
set(ax2,'Layer','top');
box on;
hold on; 

% Ramsey
tp=shadedErrorBar(rad2deg(S_ramsey.az),S_ramsey.dBdx_eq,S_ramsey.dBdx_eq_se,'k');
tp.mainLine.Color=0.3*ones(1,3);      %     'k';    % 'none';
tp.mainLine.LineWidth=1;
tp.mainLine.LineStyle='--';
tp.edge(1).Visible='off';
tp.edge(2).Visible='off';

% gradiometry
tp=shadedErrorBar(rad2deg(az),dBdx_eq,dBdxerr_eq,'r');
tp.mainLine.Color=ctheme(idx_col,:);    % 'none';
tp.mainLine.LineWidth=1;
tp.patch.FaceColor=ctheme(idx_col,:);       %cltheme(idx_col,:);
tp.patch.FaceAlpha=0.33;
tp.edge(1).Visible='off';
tp.edge(2).Visible='off';

% annotation
% zoom to theta=pi/2
[~,iaz_pi2]=min(abs(az-pi/2));
dBdx_pi2=dBdx_eq(iaz_pi2);
xlim(90+10*[-1,1]);
ylim(dBdx_pi2*(1+1*[-1,1]));
% set(ax2,'XTickLabel',[]);
ax2.FontSize=ax.FontSize-1;


%% save output
%TODO
% vars_to_save={'az','el','iel_0','gaz','gel',...
%     'tau','P_k_avg','P_k_se',...
%     'dBdx_eq_int','dBdx_eq_se_int',...
%     'B_Pramsey_halo','Berr_Pramsey_halo',...
%     'dBdx','Berrk_Pramsey',...
%     'Bk_eq','Bkerr_eq',...
%     'Br','Brerr'};  
% 
% save(['out_',getdatetimestr,'.mat'],vars_to_save{:});



%% End of script
toc