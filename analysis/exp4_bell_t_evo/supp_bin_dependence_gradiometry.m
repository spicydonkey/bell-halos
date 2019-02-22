%% Time evolution of Psi+ halo
% DKS
% 2018-10-25

tic

%% CONFIGS
% LOAD ------------------------------------------------
%   - normalised momentum vectors
%   - categorised in t-delay
configs.fdata='C:\Users\HE BEC\Dropbox\PhD\projects\halo_metrology\analysis\exp4_tevo\exp4_20181029.mat';


% EXP SCHEME  ------------------------------------------
configs.exp.v_sep=120e-3;           % pair separation velocity [m/s]
configs.exp.r_bec=50e-6;            % BEC TF radius [m]

% ramsey 
configs.exp.t_ramsey=3e-3;          % ramsey separation time [s]
configs.exp.d_sep_ramsey=configs.exp.t_ramsey * configs.exp.v_sep;  % pair separation [m]

% gradiometry
configs.exp.t_grad=1.7e-3;          % gradiometry separation time [s] (exp: 0-1.7 ms)
configs.exp.d_sep_grad=configs.exp.t_grad * configs.exp.v_sep;

% far-field spatial uncertainty (dist/angular units)
configs.exp.sig_r=configs.exp.r_bec/sqrt(2);   
configs_exp.sig_beta_ramsey=configs.exp.sig_r/(configs.exp.d_sep_ramsey/2);  
configs.exp.sig_beta_grad=configs.exp.r_bec/(configs.exp.d_sep_grad/2);    


% GRIDS/BINS ------------------------------------------------------
% grids around sphere
configs.bins.az_lim=[0,pi];
configs.bins.el_lim=[0,0];      % pi/4*[-1,1]

configs.bins.n_az=2;
configs.bins.n_el=1;

% half-cone angle
configs.bins.alpha_scan=pi*logspace(-1,0.5,20);     % min: atom num limited; max: asymptote at pi/2 


% g2 -------------------------------------------------------------
configs.g2.n_dk=7;
configs.lim_dk=[-0.2,0.2];


% bootstrapping --------------------------------------------------
configs.bootstrap.samp_frac=0.2;
configs.bootstrap.n_rep=10;


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

config_fig.idx_col=2;   % main color

% misc
f_units='normalized';
f_pos=[0.2,0.2,0.2,0.3];
f_pos_wide=[0.2,0.2,0.25,0.3];

mark_siz=7;
line_wid=1.5;
fontsize=12;
ax_lwidth=1.2;


%% load preprocessed data
load(configs.fdata);


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
az=linspace(configs.bins.az_lim(1),configs.bins.az_lim(2),configs.bins.n_az+1);
az=az(1:end-1);       % exclude last
el=linspace(configs.bins.el_lim(1),configs.bins.el_lim(2),configs.bins.n_el);
% el=0;          % THIS LINE IS DEV investigation


[gaz,gel]=ndgrid(az,el);    % AZ-EL grid
n_zone=numel(gaz);

% % special angles
% [~,iaz_disp]=arrayfun(@(x) min(abs(az-x)),az_disp);     % idx to ~displayable azim angles
% [~,iel_0]=min(abs(el));        % idx to ~zero elev angle (equator)


%% create g2 params
dk_ed_vec=linspace(configs.lim_dk(1),configs.lim_dk(2),configs.g2.n_dk+1);
dk_cent_vec=dk_ed_vec(1:end-1)+0.5*diff(dk_ed_vec);
[~,idx_dk0]=min(abs(dk_cent_vec));    % idx bin nearest zero
dk_ed={dk_ed_vec,dk_ed_vec,dk_ed_vec};
dk_grid_size=cellfun(@(k) length(k)-1,dk_ed);


%% LOOP analysis
% initialise
temp_nan=NaN(configs.bins.n_az,configs.bins.n_el,length(configs.bins.alpha_scan));
temp_nan2=NaN(length(k_tau),configs.bins.n_az,configs.bins.n_el,length(configs.bins.alpha_scan));
% temp_cell=cell(length(configs.bins.alpha_scan),1);

v_beta=temp_nan;
v_beta_se=temp_nan;

v_dBdx=temp_nan;
v_dBdx_se=temp_nan;

% v_fit={};
% v_Pi0_fit={};

v_Pi0 = temp_nan2;
v_Pi0_bs_se = temp_nan2;


H_parity=figure('Name','triplet_halo_tevo','Units',f_units,'Position',f_pos,'Renderer',config_fig.rend);

% number crunching!
for i_alpha = 1:numel(configs.bins.alpha_scan)
    alpha = configs.bins.alpha_scan(i_alpha);
    %% g2 ===============================================================
    n_tau=length(k_tau);
    
    % preallocate
    c4=cell(3,n_tau,configs.bins.n_az,configs.bins.n_el);     % DIM: corr X tau X theta X phi
    c3=cell(n_tau,configs.bins.n_az,configs.bins.n_el);       % DIM: tau X theta X phi
    m4=NaN(3,n_tau,configs.bins.n_az,configs.bins.n_el);
    m3=NaN(n_tau,configs.bins.n_az,configs.bins.n_el);
    
    g2=c4;
    g0=m4;
    Pi=m3;
    Pi0=m3;
    
    g2_bs_mu=c4;
    g2_bs_se=c4;
    g0_bs_mu=m4;
    g0_bs_se=m4;
    Pi_bs_mu=m3;
    Pi0_bs_mu=m3;
    Pi_bs_se=m3;
    Pi0_bs_se=m3;
    
    for ii=1:n_tau
        %%
        k=k_tau{ii};
        n_data_size=size(k,1);
        
        progressbar(0);
        for jj=1:n_zone
            %%
            [iaz,iel]=ind2sub([configs.bins.n_az,configs.bins.n_el],jj);
            taz=gaz(jj);
            tel=gel(jj);
            
            % capture vecs in section
            %         [~,b_A]=cellfun(@(x) inCone(x,taz,tel,alpha),k,'UniformOutput',false);     % in k
            %         [~,b_B]=cellfun(@(x) inCone(x,taz+pi,-tel,alpha),k,'UniformOutput',false);    % in -k
            %         b_AB=cellfun(@(b1,b2) b1|b2,b_A,b_B,'UniformOutput',false);        % in both
            %         k_in=cellfun(@(x,b) x(b,:),k,b_AB,'UniformOutput',false);
            k_in=cellfun(@(x) inDoubleCone(x,taz,tel,alpha),k,'UniformOutput',false);
            
            %% g2
            tg2=summary_disthalo_g2(k_in,dk_ed,0,0,0,0);
            
            tg0=cellfun(@(g) g(idx_dk0,idx_dk0,idx_dk0),tg2);    % BB corr strength
            [tPi,tPi0]=g2toE(mean(tg0(1:2)),tg0(3));      % correlator
            
            %%% store
            g2(:,ii,iaz,iel)=tg2(:);
            g0(:,ii,iaz,iel)=tg0;
            Pi(ii,iaz,iel)=tPi;
            Pi0(ii,iaz,iel)=tPi0;
            
            %% BOOTSTRAPPING
            % set up
            bs_nsamp=round(configs.bootstrap.samp_frac*n_data_size);    % # data to sample for bs
            bs_Isamp=cellfun(@(c) randi(n_data_size,[bs_nsamp,1]), cell(configs.bootstrap.n_rep,1),...
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
            tg2_bs_se=cellfun(@(x) sqrt(configs.bootstrap.samp_frac)*std(x,0,4),tg2_bs_cat,'UniformOutput',false);
            
            tg0_bs_mu=mean(tg0_bs,1);
            tg0_bs_se=sqrt(configs.bootstrap.samp_frac)*std(tg0_bs,0,1);
            
            tB_bs_mu=mean(tB_bs,1);
            tB0_bs_mu=mean(tB0_bs,1);
            tB_bs_se=sqrt(configs.bootstrap.samp_frac)*std(tB_bs,0,1);
            tB0_bs_se=sqrt(configs.bootstrap.samp_frac)*std(tB0_bs,0,1);
            
            %%% store
            g2_bs_mu(:,ii,iaz,iel)=tg2_bs_mu(:);
            g2_bs_se(:,ii,iaz,iel)=tg2_bs_se(:);
            g0_bs_mu(:,ii,iaz,iel)=tg0_bs_mu;
            g0_bs_se(:,ii,iaz,iel)=tg0_bs_se;
            Pi_bs_mu(ii,iaz,iel)=tB_bs_mu;
            Pi0_bs_mu(ii,iaz,iel)=tB0_bs_mu;
            Pi_bs_se(ii,iaz,iel)=tB_bs_se;
            Pi0_bs_se(ii,iaz,iel)=tB0_bs_se;
            
            progressbar(jj/n_zone);
        end
    end
    
    %% store result -----------------------------------------------
    v_Pi0(:,:,:,i_alpha) = Pi0;
    v_Pi0_bs_se(:,:,:,i_alpha) = Pi0_bs_se;
    
    
    %% vis: Parity - ALL
    %%% ALL
    [cc,ccl,ccd]=palette(n_zone);

    figure(H_parity);
    hold on;
    
    for ii=1:n_zone
        [iaz,iel]=ind2sub([configs.bins.n_az,configs.bins.n_el],ii);
        
        tp=ploterr(tau,Pi0(:,iaz,iel),[],Pi0_bs_se(:,iaz,iel),'-o','hhxy',0);
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
    %      B(x) = Pi0 + dBdx * x    --- const grad mag-field
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
    mdl_tevo2.fit=cell(configs.bins.n_az,configs.bins.n_el);
    for ii=1:n_zone
        [iaz,iel]=ind2sub([configs.bins.n_az,configs.bins.n_el],ii);
        
        tPi0=Pi0(:,iaz,iel);
        mdl_tevo2.fit{iaz,iel}=fitnlm(tau,tPi0,mdl_tevo2.mdl,mdl_tevo2.par0,...
            'CoefficientNames',mdl_tevo2.cname,'Options',mdl_tevo2.fopt);
    end
    
    % eval fitted model
    t_fit=linspace(0,5,1e3);     % time since collision pulse to eval fit
    mdl_tevo2.fit_par=arrayfun(@(I) cellfun(@(m) m.Coefficients.Estimate(I),mdl_tevo2.fit),1:numel(mdl_tevo2.cname),'UniformOutput',false);
    mdl_tevo2.fit_par_se=arrayfun(@(I) cellfun(@(m) m.Coefficients.SE(I),mdl_tevo2.fit),1:numel(mdl_tevo2.cname),'UniformOutput',false);
    Pi0_fit2=cellfun(@(f) feval(f,t_fit),mdl_tevo2.fit,'UniformOutput',false);
    
    
    %%% Magnetic field gradient: dBdx -----------------------------------------
    beta=abs(mdl_tevo2.fit_par{1});             % get fit params and ensure sign is positive
    beta_se=abs(mdl_tevo2.fit_par_se{1});
    
    
    % far-field
    dBdx=1e6*beta/(configs.exp.v_sep*C_gymag/2);            % [G/m] (factor to convert time-scaling ms^-2 --> s^-2)
    dBdx_se=1e6*beta_se/(configs.exp.v_sep*C_gymag/2);
    
    
    %% STAT: store result =========================================
    v_beta(:,:,i_alpha)=beta;
    v_beta_se(:,:,i_alpha)=beta_se;
    
    v_dBdx(:,:,i_alpha)=dBdx;
    v_dBdx_se(:,:,i_alpha)=dBdx_se;
    
%     v_fit{end+1} = mdl_tevo2;
%     v_Pi0_fit{end+1} = Pi0_fit2;
end

%% DATA VIS ============================================================
%% VIS: parity vs bin-size
h=figure('Name','binsize_vs_parity','Units',config_fig.units,'Position',[0,0,8.6,4.5*3],'Renderer',config_fig.rend);

tau_plot=round(linspace(1,n_tau,3));

for kk=1:length(tau_plot)
    itau=tau_plot(kk);
    subplot(3,1,kk);
    ax=gca;
    hold on;
    
    pleg=[];
    for ii=1:configs.bins.n_az
        for jj=1:configs.bins.n_el
            tI=sub2ind([configs.bins.n_az,configs.bins.n_el],ii,jj);
            
            tp=shadedErrorBar(configs.bins.alpha_scan/pi,squeeze(v_Pi0(itau,ii,jj,:)),squeeze(v_Pi0_bs_se(itau,ii,jj,:)));
            tp.mainLine.LineStyle=config_fig.line_sty{ii};
            tp.mainLine.DisplayName=sprintf('%0.3g, %0.3g',gaz(ii),gel(jj));
            tp.patch.FaceAlpha=0.33;
            tp.edge(1).LineStyle=config_fig.line_sty{ii};
            tp.edge(2).LineStyle=config_fig.line_sty{ii};
            
            pleg(end+1)=tp.mainLine;
        end
    end

    
    title(sprintf('%0.2g ms',tau(itau)));
    
    % annotation
    set(ax,'Layer','top');
    set(ax,'XScale','log');
    xlabel('bin size $\alpha/\pi$');
    ylabel('Parity');
    box on;
    ax.FontSize=config_fig.ax_fontsize;
    ax.LineWidth=config_fig.ax_lwid;
    
    ax.XLim(1)=0.1;     % below this not enough data for g2
    ylim(1.5*[-1, 1]);
    
    
    % bin cut-off ----------------------
    alpha_max=0.5*pi;        % cutoff for entanglement analysis
    xdata_cutoff=[alpha_max/pi,ax.XLim(2),ax.XLim(2),alpha_max/pi];
    col_cutoff = 0.9*ones(1,3);
    
    o_cutoff=patch('XData',xdata_cutoff,'YData',[ax.YLim(1),ax.YLim(1),ax.YLim(2),ax.YLim(2)],...
        'FaceColor',col_cutoff,'EdgeColor','none');
    uistack(o_cutoff,'bottom');
    
    if kk==1
        lgd=legend(pleg);
        lgd.Title.String='$\theta, \phi$';
        lgd.Location='southeast';
        lgd.Box='off';
    end
end



%% VIS: ERR parity vs bin-size
h=figure('Name','binsize_vs_parity_err','Units',config_fig.units,'Position',[0,0,8.6,4.5*3],'Renderer',config_fig.rend);

tau_plot=round(linspace(1,n_tau,3));

for kk=1:length(tau_plot)
    itau=tau_plot(kk);
    subplot(3,1,kk);
    ax=gca;
    hold on;
    
    pleg=[];
    for ii=1:configs.bins.n_az
        for jj=1:configs.bins.n_el            
            tI=sub2ind([configs.bins.n_az,configs.bins.n_el],ii,jj);
            
            tp=plot(configs.bins.alpha_scan/pi,squeeze(v_Pi0_bs_se(itau,ii,jj,:)));
            tp.LineStyle=config_fig.line_sty{ii};
            tp.Marker='none';
            tp.Color='k';
            tp.DisplayName=sprintf('%0.3g, %0.3g',gaz(ii),gel(jj));
            
            
            pleg(end+1)=tp;
        end
    end

    
    title(sprintf('%0.2g ms',tau(itau)));
    
    % annotation
    set(ax,'Layer','top');
    set(ax,'XScale','log');
    set(ax,'YScale','log');
    xlabel('bin size $\alpha/\pi$');
    ylabel('Unc Parity');
    box on;
    ax.FontSize=config_fig.ax_fontsize;
    ax.LineWidth=config_fig.ax_lwid;
    
    ax.XLim(1)=0.1;     % below this not enough data for g2
%     ylim(1.5*[-1, 1]);
    
    % bin cut-off ----------------------
    alpha_max=0.5*pi;        % cutoff for entanglement analysis
    xdata_cutoff=[alpha_max/pi,ax.XLim(2),ax.XLim(2),alpha_max/pi];
    col_cutoff = 0.9*ones(1,3);
    
    o_cutoff=patch('XData',xdata_cutoff,'YData',[ax.YLim(1),ax.YLim(1),ax.YLim(2),ax.YLim(2)],...
        'FaceColor',col_cutoff,'EdgeColor','none');
    uistack(o_cutoff,'bottom');
    
    % legend---------------------------
    if kk==1
        lgd=legend(pleg);
        lgd.Title.String='$\theta, \phi$';
        lgd.Location='northeast';
        lgd.Box='off';
    end
end


% TODO: if run out of linestyles, do with colors?
%
%% VIS: binsize vs mean
h=figure('Name','binsize_vs_dBdx','Units',config_fig.units,'Position',[0,0,8.6,4.5],'Renderer',config_fig.rend);
ax=gca;
hold on;

pleg=[];
for ii=1:configs.bins.n_az
    for jj=1:configs.bins.n_el
        tI=sub2ind([configs.bins.n_az,configs.bins.n_el],ii,jj);

        tp=shadedErrorBar(configs.bins.alpha_scan/pi,squeeze(v_dBdx(ii,jj,:)),squeeze(v_dBdx_se(ii,jj,:)));
        tp.mainLine.LineStyle=config_fig.line_sty{ii};
        tp.mainLine.DisplayName=sprintf('%0.3g, %0.3g',gaz(ii),gel(jj));
        tp.patch.FaceAlpha=0.33;

        pleg(end+1)=tp.mainLine;
    end
end

% annotation
set(ax,'Layer','top');
set(ax,'XScale','log');
xlabel('bin size $\alpha/\pi$');
ylabel('$dB/dx$ (G/m)');
box on;
ax.FontSize=config_fig.ax_fontsize;
ax.LineWidth=config_fig.ax_lwid;

ax.XLim(1)=0.1;     % below this not enough data for g2
ylim([0, 5]);


% bin cut-off ----------------------
alpha_max=0.5*pi;        % cutoff for entanglement analysis
xdata_cutoff=[alpha_max/pi,ax.XLim(2),ax.XLim(2),alpha_max/pi];
col_cutoff = 0.9*ones(1,3);

o_cutoff=patch('XData',xdata_cutoff,'YData',[ax.YLim(1),ax.YLim(1),ax.YLim(2),ax.YLim(2)],...
    'FaceColor',col_cutoff,'EdgeColor','none');
uistack(o_cutoff,'bottom');

% legend ----------------------
lgd=legend(pleg);
lgd.Title.String='$\theta, \phi$';
lgd.Location='northeast';
lgd.Box='off';



%% VIS: binsize vs error
h=figure('Name','binsize_vs_dBdx_err','Units',config_fig.units,'Position',[0,0,8.6,4.5],'Renderer',config_fig.rend);
ax=gca;
hold on;
for ii=1:configs.bins.n_az
    for jj=1:configs.bins.n_el
        tI=sub2ind([configs.bins.n_az,configs.bins.n_el],ii,jj);
        tp=plot(configs.bins.alpha_scan/pi,squeeze(v_dBdx_se(ii,jj,:)));
%         tp=plot(configs.bins.alpha_scan,squeeze(Berr_alpha(ii,jj,:)));
        tp.LineStyle=config_fig.line_sty{ii};
        tp.Marker='none';
        tp.Color='k';
    end
end

% annotation
set(ax,'Layer','top');
set(ax,'XScale','log');
set(ax,'YScale','log');
xlabel('bin size $\alpha/\pi$');
ylabel('Unc $dB/dx$ (G/m)');
box on;
ax.FontSize=config_fig.ax_fontsize;
ax.LineWidth=config_fig.ax_lwid;

ax.XLim(1)=0.1;     % below this not enough data for g2

% bin cut-off ----------------------
alpha_max=0.5*pi;        % cutoff for entanglement analysis
xdata_cutoff=[alpha_max/pi,ax.XLim(2),ax.XLim(2),alpha_max/pi];
col_cutoff = 0.9*ones(1,3);

o_cutoff=patch('XData',xdata_cutoff,'YData',[ax.YLim(1),ax.YLim(1),ax.YLim(2),ax.YLim(2)],...
    'FaceColor',col_cutoff,'EdgeColor','none');
uistack(o_cutoff,'bottom');


%% save
vars_to_save={'configs','config_fig',...
    'v_beta','v_beta_se','v_dBdx','v_dBdx_se',...
    'v_Pi0', 'v_Pi0_bs_se',...
    'tau',...
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

