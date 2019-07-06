% Parity analysis
% phase uncertainty analysis
%
% 2019-07-05
% DKS
%


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
% configs.bins.az_lim=[0,pi];
% configs.bins.el_lim=pi/4*[-1,1];      
% 
% configs.bins.n_az=60;
% configs.bins.n_el=30;

% bin size (half-cone angle)
configs.bins.alpha=configs.exp.sig_beta_grad;      

% misc
% configs.bins.az_disp=deg2rad([0,45,90,135]);

% vis ----------------------------------------------------------------
% publication
config_fig = loadFigureConfig;

config_fig.mark_typ={'+','o','^','d'};
config_fig.line_sty={'-','--',':','-.'};

config_fig.col_theme=parula(5);       % theme color
config_fig.coll_theme=colshades(config_fig.col_theme);


%% load preprocessed data
load(configs.fdata);


%% transform to RH coords
% SAVE ORIGINAL 
if ~exist('k_tau_orig','var')
    k_tau_orig=k_tau;
end

k_tau=cellfun(@(C) cellfun(@(x) tzxy2RHtzxy2(x),C,'UniformOutput',false),...
    k_tau_orig,'UniformOutput',false);      % EXP-coord sys (z against g)


clearvars('k_tau_orig');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% MAIN
% prototyping stage
%

% preprocess input data
idx_tau = 5;        % 5 is 1.4ms and low mode occupancy
k = k_tau{idx_tau};        


% config analysis
alpha_intg = configs.bins.alpha;         % integration volume
alpha_corr = 0.06;                       % correlator volume

el_scan = 0;
az_scan = 0:alpha_intg:pi;


%%% unit analysis - for specific az,el angle
el_intg = el_scan;
% az_intg = az_scan(end);
az_intg = 0;


% config corr volumes
daz_corr = 1.5*alpha_corr;
del_corr = 1.5*alpha_corr;

% simple azel-rectangular grid search: problems at boundary and rectangular
%   grid search at corners
% az_corr_scan = (az_intg-alpha_intg):daz_corr:(az_intg+alpha_intg);
% el_corr_scan = (el_intg-alpha_intg):del_corr:(el_intg+alpha_intg);

% simple azel-rectangular scan centered with intg vol
az_corr_scan = 0:daz_corr:alpha_intg;
az_corr_scan = az_intg + [-az_corr_scan(1:end-1),az_corr_scan];
el_corr_scan = 0:del_corr:alpha_intg;
el_corr_scan = el_intg + [-el_corr_scan(1:end-1),el_corr_scan];

naz_corr_scan = length(az_corr_scan);
nel_corr_scan = length(el_corr_scan);

[gaz_corr_scan, gel_corr_scan] = ndgrid(az_corr_scan,el_corr_scan);     %grid
ngrid = numel(gaz_corr_scan);

% reduce to atoms in integration vol
tk_intg = cellfun(@(x) inDoubleCone(x,az_intg,el_intg,alpha_intg),k,'UniformOutput',false);

N_intg = shotSize(tk_intg);

% search corr volumes
parity = cell(naz_corr_scan,nel_corr_scan);

for ii=1:naz_corr_scan
    taz_corr = az_corr_scan(ii);
    for jj=1:nel_corr_scan    
        tel_corr = el_corr_scan(jj);

        if diffAngleSph(az_intg,el_intg,taz_corr,tel_corr) <= alpha_intg - alpha_corr    
            % evaluate correlator/parity
            n_A=cellfun(@(x) size(inCone(x,taz_corr,tel_corr,alpha_corr),1),tk_intg);
            n_B=cellfun(@(x) size(inCone(x,taz_corr+pi,-tel_corr,alpha_corr),1),tk_intg);
            
            S_A = diff(n_A,1,2);
            S_B = diff(n_B,1,2);
            
            tparity = S_A .* S_B;
            
            N_corr = n_A + n_B;
            dN = N_intg - N_corr;
            
            % post-selection for (NA = 1, NB = 1) pair detection events
            tparity_postsel = tparity;
            tparity_postsel(~(abs(tparity)==1)) = NaN;
        else
            tparity=[];
        end
        
        parity{ii,jj} = tparity;
    end
end

parity_collated = cat(1,parity{:});

parity_collated_ps = parity_collated;
parity_collated_ps(~(abs(parity_collated)==1)) = NaN;


parity_mean = mean(parity_collated_ps,'omitnan');
parity_sdev = std(parity_collated_ps,'omitnan');

phase_mean = 0.5*acos(parity_mean);
phase_sdev = parity_sdev/(2*sin(2*phase_mean));

%%% bootstrapping
% to estimate "statistical uncertainty" in:
n_data=numel(parity_collated_ps);
B = 200;        % bootstrap repetition
% 1. mean parity
idx_bssamp = randi(n_data,B,n_data);
parity_bssamp = parity_collated_ps(idx_bssamp);     % bootstrap samples
parity_mean_bs = mean(parity_bssamp,2,'omitnan');   % mean from bs
parity_mean_bserr = std(parity_mean_bs);            % evaluate uncertainty

% 2. sdev parity
idx_bssamp = randi(n_data,B,n_data);
parity_bssamp = parity_collated_ps(idx_bssamp);     % bootstrap samples
parity_sdev_bs = std(parity_bssamp,0,2,'omitnan');  
parity_sdev_bserr = std(parity_sdev_bs);            % evaluate uncertainty

% 3. mean phase
idx_bssamp = randi(n_data,B,n_data);
parity_bssamp = parity_collated_ps(idx_bssamp);     % bootstrap samples
parity_mean_bs = mean(parity_bssamp,2,'omitnan');   % mean from bs
phase_mean_bs = 0.5*acos(parity_mean_bs);           % phase from parity mean
phase_mean_bserr = std(phase_mean_bs);           % evaluate uncertainty

% 3. sdev phase
idx_bssamp = randi(n_data,B,n_data);
parity_bssamp = parity_collated_ps(idx_bssamp);     % bootstrap samples
% parity_mean_bs = mean(parity_bssamp,2,'omitnan');   
parity_sdev_bs = std(parity_bssamp,0,2,'omitnan');  
% phase_mean_bs = 0.5*acos(parity_mean_bs);           % phase from parity mean
% phase_sdev_bs = parity_sdev_bs./(2*sin(2*phase_mean_bs));
phase_sdev_bs = parity_sdev_bs./(2*sin(2*phase_mean));  % parity sensitivity from phi estimation
phase_sdev_bserr = std(phase_sdev_bs);


% summarise bootstrapped unc as fractional unc
parity_mean_ferr = 1e2* abs(parity_mean_bserr/parity_mean);
parity_sdev_ferr = 1e2* abs(parity_sdev_bserr/parity_sdev);
phase_mean_ferr = 1e2* abs(phase_mean_bserr/phase_mean);
phase_sdev_ferr = 1e2* abs(phase_sdev_bserr/phase_sdev);


% % phase estimation by single-shot
% % NOTE: does not seem legitimate based on the framework of quantum metrology
% phase_single = 0.5*acos(parity_collated_ps);
% phase_single_mean = mean(phase_single,'omitnan');
% phase_single_sdev = std(phase_single,'omitnan');


%% plot
H=figure('Name','embg_phi_est');
H.Renderer='painters';

JJ_max = max(abs(parity_collated));

%%% all coincidence events (JJ)
subplot(2,1,1);

JJ_edge = -JJ_max-0.5:1:JJ_max+0.5;
JJ_cent = edge2cent(JJ_edge);
N_JJ = histcounts(parity_collated,JJ_edge);
Ntot_JJ = sum(N_JJ);

P_JJ = N_JJ/Ntot_JJ;        % probability
P_JJ_sdev_est = sqrt(P_JJ.*(1-P_JJ)/Ntot_JJ);   % err estimate

% barplot
p_bar_JJ = bar(JJ_cent,P_JJ,'b');            
p_bar_JJ.FaceAlpha = 0.5;

hold on;

% errorbar
% TODO: check calculation
p_er_JJ = errorbar(JJ_cent,P_JJ,P_JJ_sdev_est,P_JJ_sdev_est);    
p_er_JJ.Color = [0 0 0];                            
p_er_JJ.LineStyle = 'none';  

titlestr = sprintf('$\\tau=%0.2f$ ms; $(\\theta,\\phi)=(%0.2g, %0.2g)$ deg; $(\\alpha,\\delta k) = (%0.2f,\\,%0.2f)$; $N=%d$',tau(idx_tau),rad2deg(az_intg),rad2deg(el_intg),alpha_intg,alpha_corr,Ntot_JJ);
title(titlestr);
xlabel('$J_x^A J_x^B$');
ylabel('probability');

set(gca,'YScale','log');
set(gca,'YMinorTick','off');
set(gca,'XTick',-5:5);

%%% post-selected
subplot(2,1,2);

% histogram(parity_collated_ps);

ss_edge = -1.5:1:1.5;
ss_cent = edge2cent(ss_edge);
N_ss = histcounts(parity_collated_ps,ss_edge);
Ntot_ss = sum(N_ss);

P_ss = N_ss/Ntot_ss;        % probability
P_ss_sdev_est = sqrt(P_ss.*(1-P_ss)/Ntot_ss);   % err estimate

% barplot
p_bar_ss = bar(ss_cent,P_ss,'r');            
p_bar_ss.FaceAlpha = 0.5;

hold on;

% errbar
% TODO: estimate err
p_er_ss = errorbar(ss_cent,P_ss,P_ss_sdev_est,P_ss_sdev_est);    
p_er_ss.Color = [0 0 0];                            
p_er_ss.LineStyle = 'none';  


set(gca,'XTick',-5:5);
titlestr = sprintf('post-selected for single pair: $N=%d$',Ntot_ss);
title(titlestr);
xlabel('parity $\sigma_x^A \sigma_x^B$')
ylabel('probability');


fig_text_parity = sprintf('$\\textbf{parity}$\navg = %0.2f $\\pm$ %0.1g\ns.d. = %0.2f $\\pm$ %0.1g',...
    parity_mean,round(parity_mean_bserr,1,'significant'),parity_sdev,round(parity_sdev_bserr,1,'significant'));
fig_text_phase = sprintf('\n$\\textbf{phase}\\,\\Phi$\navg = %0.2f $\\pm$ %0.1g\ns.d. = %0.2f $\\pm$ %0.1g',...
    phase_mean,round(phase_mean_bserr,1,'significant'),phase_sdev,round(phase_sdev_bserr,1,'significant'));
fig_text=[fig_text_parity,fig_text_phase];
text(0,mean(get(gca,'YLim')),fig_text ,'HorizontalAlignment','center');


