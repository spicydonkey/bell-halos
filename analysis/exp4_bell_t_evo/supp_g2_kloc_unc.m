%% SOM ANALYSIS: uncertainty in localised g2 by bootstrapping
% DKS
% 2019-03-27


%% CONFIGS -------------------------------------------------
% data file
%   contains spin,momentum resolved vectors
configs.data.main='C:\Users\HE BEC\Dropbox\PhD\projects\halo_metrology\analysis\exp4_tevo\exp4_20181029.mat';
configs.data.supp = 'C:\Users\HE BEC\Documents\MATLAB\bell-halos\analysis\exp4_bell_t_evo\out\supp_kloc_g2_vs_tau_20190326_hr.mat';

% bootstrapping --------------------------------------------------
configs.bootstrap.B=20;

% g2
% single-bin at dk=0
configs.g2.dk_lim=0.05*[-1,1];
configs.g2.dk_n=1;

% k-modes
configs.mode.alpha=pi/10;           % cone half-angle

% k-modes around equator
configs.mode.n_azel=11;
configs.mode.azel=linspace(-pi,pi,configs.mode.n_azel)';        % azim (rad)
configs.mode.azel(:,2)=0*configs.mode.azel(:,1);                % elev (rad)

% vis
config_fig = loadFigureConfig;      % load template


%% load data -------------------------------------------------
load(configs.data.main);

% SUPP data: load high-res g2 profiles
S_data = load(configs.data.supp);


%% debug: reduce data
% warning('Reducing data for speed.');
% k_tau{1}=k_tau{1}(1:3000,:);
% for ii=1:numel(k_tau)
%     k_tau{ii}=k_tau{ii}(1:round(length(k_tau{ii})/sqrt(30)),:);
% end

% cull tau
% warning('culling Tau');
% k_tau=k_tau(6);


%% MAIN ------------------------------------------------------
% % configure g2 analysis
dk_ed_vec=linspace(configs.g2.dk_lim(1),configs.g2.dk_lim(2),configs.g2.dk_n+1);
dk_cent_vec=edge2cent(dk_ed_vec);
% [~,idx_dk0]=min(abs(dk_cent_vec));    % idx bin nearest zero
dk_ed={dk_ed_vec,dk_ed_vec,dk_ed_vec};
dk_grid_size=cellfun(@(k) length(k)-1,dk_ed);


% % configure k-mode
n_mode=size(configs.mode.azel,1);


% % localised g2 analysis
n_tau=length(k_tau);
n_shots=shotSize(k_tau);

% preallocate
g2_mode_bs=cell(n_tau,1);

g2_mode_avg = cell(n_tau,1);
g2_mode_se = cell(n_tau,1);

for ii=1:n_tau
    % preallocate
    g2_mode_avg{ii}=NaN(n_mode,3);
    g2_mode_se{ii}=NaN(n_mode,3);
    
    g2_mode_bs=cell(n_mode,1);
    
    t_k = k_tau{ii};
    
    % g2 spatial distribution
    progressbar(0)
    for jj=1:n_mode
        % get this mode info
        taz=configs.mode.azel(jj,1);
        tel=configs.mode.azel(jj,2);
        
        % filter atoms in this mode
        tk_mode = cellfun(@(x) inDoubleCone(x,taz,tel,configs.mode.alpha),t_k,'uni',0);
        
        % % g2 bootstrapping
        % set up bootstrapping samples
        bs_Isamp=cellfun(@(c) randi(n_shots(ii),[n_shots(ii),1]), cell(configs.bootstrap.B,1),...
            'uni',0);
        
        % run analysis
        tg2_bs=cellfun(@(I) summary_disthalo_g2(tk_mode(I,:),dk_ed,0,0,0,0),bs_Isamp,...
            'uni',0);
        
        % tidy variable structure
        tg2_bs=cell2mat(cat(1,tg2_bs{:}));    % BS-idx x g2 type: g2(0)
        
        % evaluate uncertainty
        tg2_bs_avg=mean(tg2_bs,1);
        tg2_bs_se=std(tg2_bs,1);
        
        % store
        g2_mode_avg{ii}(jj,:)=tg2_bs_avg;
        g2_mode_se{ii}(jj,:)=tg2_bs_se;
        
        g2_mode_bs{ii}{jj} = tg2_bs;	% bootstrapped sample analysis results
        
        
        progressbar(jj/n_mode);
    end
end

%% post processing
% simplify g2 corr types
g2_corr_avg=cellfun(@(x) mean(x(:,1:2),2),g2_mode_avg,'uni',0);
g2_corr_se=cellfun(@(x) vnorm(x(:,1:2),2),g2_mode_se,'uni',0);

g2_anti_avg=cellfun(@(x) x(:,3),g2_mode_avg,'uni',0);
g2_anti_se=cellfun(@(x) x(:,3),g2_mode_se,'uni',0);


%% VIS: Equatorial g2 profile
H={};
for ii=1:n_tau
   H{ii}=figure('Name','g2_kloc_equator','Units',S_data.config_fig.units,'Position',[0 0 7 3],'Renderer',S_data.config_fig.rend);
   H{ii}.Name=strcat(H{ii}.Name,'_',num2str(ii));

   hold on;
   % supps data - g2 from all data (labels reversed!)
   tp_corr = plot(rad2deg(S_data.vaz_filled),S_data.g2_anti{ii}(:,S_data.idx_el0),'r','DisplayName','$\uparrow\uparrow/\downarrow\downarrow$');
   tp_anti = plot(rad2deg(S_data.vaz_filled),S_data.g2_corr{ii}(:,S_data.idx_el0),'b--','DisplayName','$\uparrow\downarrow$');

   % bootstrapping - avg and SE
   tp_corr_bs=ploterr(rad2deg(configs.mode.azel(:,1)),g2_corr_avg{ii},[],g2_corr_se{ii},'or');
   set(tp_corr_bs(1),'Color','r','MarkerFaceColor','r',...
       'MarkerSize',config_fig.mark_siz,'LineWidth',config_fig.line_wid);
   set(tp_corr_bs(2),'Color','r','LineWidth',config_fig.line_wid);

   tp_anti_bs=ploterr(rad2deg(configs.mode.azel(:,1)),g2_anti_avg{ii},[],g2_anti_se{ii},'ob');
   set(tp_anti_bs(1),'Color','b','MarkerFaceColor','w',...
       'MarkerSize',config_fig.mark_siz,'LineWidth',config_fig.line_wid);
   set(tp_anti_bs(2),'Color','b','LineWidth',config_fig.line_wid);

   % annotate
   box on;
   ax=gca;
   ax.FontSize=S_data.config_fig.ax_fontsize;
   ax.LineWidth=S_data.config_fig.ax_lwid;
   xlabel('$\theta$ (deg)');
   ylabel('$g^{(2)}$');
   xlim([-180,180]);
   set(gca,'XTick',-180:90:180);

%    lgd=legend('Location','best');
%     lgd=legend([tp_corr,tp_anti],'Location','best');
    
%     print_vecrast(H{ii});
end