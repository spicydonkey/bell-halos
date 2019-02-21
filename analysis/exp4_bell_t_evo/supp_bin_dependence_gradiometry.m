%% Time evolution of Psi+ halo
% DKS
% 2018-10-25

tic

%% CONFIGS
% preprocessed data file
%   - normalised momentum vectors
%   - categorised in t-delay
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
% alpha=sig_psf_gradiometry;
% bin size
c_alpha=logspace(log10(0.05),log10(0.25),6);       % factors of alpha
v_alpha=pi*c_alpha;

lim_az=[0,pi];              % limit of azim angle (exc +pi)
phi_max=pi/4;
lim_el=[-phi_max,phi_max];

% n_az=40;                	% equispaced bins
% n_el=20;

n_az=1;                	% equispaced bins
n_el=1;


az_disp=deg2rad(0);     % azim sections (great circles) to display

% g2 -------------------------------------------------------------
n_dk=7;
lim_dk=[-0.2,0.2];

% bootstrapping
bs_frac=0.1;        % 0.1
bs_nrep=10;         % 20


% Physical constants ---------------------------------------------
Cphys=physConsts;
C_gymag=Cphys.He_gymag;         % He* gyromagnetic ratio (gamma) [Hz/G]

% vis ----------------------------------------------------------------
% publication
config_fig = loadFigureConfig;

config_fig.mark_typ={'+','o','^','d'};
config_fig.line_sty={'-','--',':','-.'};

config_fig.col_theme=parula(5);       % theme color
config_fig.coll_theme=colshades(config_fig.col_theme);

idx_col=2;     % mid-color

% misc
f_units='normalized';
f_pos=[0.2,0.2,0.2,0.3];
f_pos_wide=[0.2,0.2,0.25,0.3];
config_fig.rend='painters';

mark_siz=7;
line_wid=1.5;
fontsize=12;
ax_lwidth=1.2;


%% load preprocessed data
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
% Vel=linspace(lim_el(1),lim_el(2),n_el);
Vel=0;          % THIS LINE IS DEV investigation


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


%% LOOP analysis
% initialise
v_beta=[];
v_beta_se=[];

v_dBdx=[];
v_dBdx_se=[];

v_fit={};

v_Bfit={};


% vis
H_parity=figure('Name','triplet_halo_tevo','Units',f_units,'Position',f_pos,'Renderer',config_fig.rend);


% number crunching!
for i_alpha = 1:numel(v_alpha)
    alpha = v_alpha(i_alpha);
    %% g2 ===============================================================
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

%     figure(H_parity);
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
    
    drawnow;
    
    %% Gradiometry Analysis ================================================
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
    % t0 is fully constrained = 0
    mdl_tevo2.mdl='y~cos(2*beta*(x)^2)';
    mdl_tevo2.cname={'beta'};
    mdl_tevo2.fopt=statset('TolFun',10^-10,'TolX',10^-10,'MaxIter',10^6,'UseParallel',0);
    mdl_tevo2.par0=[0.15];
    
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
    
    
    
    %% STAT: store result =========================================
    v_beta(end+1)=beta;
    v_beta_se(end+1)=beta_se;
    
    v_dBdx(end+1)=dBdx;
    v_dBdx_se(end+1)=dBdx_se;
    
    v_fit{end+1} = mdl_tevo2;
    
    v_Bfit{end+1} = B0_fit2;
end

%% VIS: binsize vs mean ==================================================
h=figure('Name','binsize_vs_mean_dBdx','Units',config_fig.units,'Position',[0,0,8.6,4.5],'Renderer',config_fig.rend);
ax=gca;
hold on;

pleg=[];
for ii=1:n_az
    for jj=1:n_el
        tI=sub2ind([n_az,n_el],ii,jj);

        tp=shadedErrorBar(v_alpha/pi,v_dBdx,v_dBdx_se);
        tp.mainLine.LineStyle=config_fig.line_sty{tI};
        tp.mainLine.DisplayName=sprintf('%0.3g, %0.3g',vaz(ii),vel(jj));
        tp.patch.FaceAlpha=0.33;

        pleg(end+1)=tp.mainLine;
    end
end

lgd=legend(pleg);
lgd.Title.String='$\theta, \phi$';
lgd.Location='northeast';
lgd.Box='off';

% annotation
set(ax,'Layer','top');
set(ax,'XScale','log');
xlabel('bin size $\alpha/\pi$');
ylabel('Mean (G/m)');
box on;
ax.FontSize=config_fig.ax_fontsize;
ax.LineWidth=config_fig.ax_lwid;

%% VIS: binsize vs error
h=figure('Name','binsize_vs_error_dBdx','Units',config_fig.units,'Position',[0,0,8.6,4.5],'Renderer',config_fig.rend);
ax=gca;
hold on;
for ii=1:n_az
    for jj=1:n_el
        tI=sub2ind([n_az,n_el],ii,jj);
        tp=plot(v_alpha/pi,v_dBdx_se);
%         tp=plot(v_alpha,squeeze(Berr_alpha(ii,jj,:)));
        tp.LineStyle=config_fig.line_sty{tI};
        tp.Marker='none';
        tp.Color='k';
    end
end

% annotation
set(ax,'Layer','top');
set(ax,'XScale','log');
set(ax,'YScale','log');
xlabel('bin size $\alpha/\pi$');
ylabel('Uncertainty (G/m)');
box on;
ax.FontSize=config_fig.ax_fontsize;
ax.LineWidth=config_fig.ax_lwid;



%% save
vars_to_save={'v_alpha','c_alpha',...
    'config_fig',...
    'v_beta','v_beta_se','v_dBdx','v_dBdx_se',...
    'v_fit','v_Bfit',...
    };

% check exists
for ii=1:numel(vars_to_save)
    tvarname=vars_to_save{ii};
    if ~exist(tvarname,'var')
        warning('variable %s does not exist.',tvarname);
    end
end

% save
save(['supp_binsize_',getdatetimestr,'.mat'],vars_to_save{:});

