%% SOM ANALYSIS: localised g2
% DKS
% 2019-03-20

%% CONFIGS -------------------------------------------------
% data file
%   contains spin,momentum resolved vectors
fdata='C:\Users\HE BEC\Dropbox\PhD\projects\halo_metrology\analysis\exp4_tevo\exp4_20181029.mat';

% g2 
configs.g2.dk_lim=0.02*[-1,1];       
configs.g2.dk_n=1;                 

% k-modes
configs.mode.alpha=pi/10;                 % cone half-angle
configs.mode.lim_az=pi*[-1,1];              % limit of azim angle (exc +pi)
configs.mode.lim_el=pi/2*[-1,1];
configs.mode.n_az=20;                	% equispaced bins
configs.mode.n_el=10;

% vis 
config_fig = loadFigureConfig;      % load template

%% load data -------------------------------------------------
load(fdata);

%%% debug: reduce data
% warning('Reducing data for speed.');
% for ii=1:numel(k_tau)
%     k_tau{ii}=k_tau{ii}(1:round(length(k_tau{ii})/10),:);
% end

%% MAIN ------------------------------------------------------
% % configure g2 analysis
dk_ed_vec=linspace(configs.g2.dk_lim(1),configs.g2.dk_lim(2),configs.g2.dk_n+1);
dk_cent_vec=dk_ed_vec(1:end-1)+0.5*diff(dk_ed_vec);
[~,idx_dk0]=min(abs(dk_cent_vec));    % idx bin nearest zero
dk_ed={dk_ed_vec,dk_ed_vec,dk_ed_vec};
dk_grid_size=cellfun(@(k) length(k)-1,dk_ed);

% % configure k-modes
vaz=linspace(configs.mode.lim_az(1),configs.mode.lim_az(2),configs.mode.n_az+1);
vaz=vaz(1:end-1);       % exclude last
vel=linspace(configs.mode.lim_el(1),configs.mode.lim_el(2),configs.mode.n_el);

[gaz,gel]=ndgrid(vaz,vel);    % AZ-EL grid
n_zone=numel(gaz);

% % localised g2 analysis
n_tau=length(k_tau);

g2_mode = cell(n_tau,1);

for ii=1:n_tau
    % preallocate
    g2_mode{ii}=NaN(configs.mode.n_az,configs.mode.n_el,3);

    t_k = k_tau{ii};    

    % g2 spatial distribution
    progressbar(0)
    for jj=1:n_zone
        [iaz,iel]=ind2sub(size(gaz),jj);
        taz=gaz(jj);
        tel=gel(jj);

        tk_mode = cellfun(@(x) inDoubleCone(x,taz,tel,configs.mode.alpha),t_k,'uni',0);

        % g2
        t_g2 = summary_disthalo_g2(tk_mode,dk_ed,0,0,0,0);

        g2_mode{ii}(iaz,iel,:) = [t_g2{:}];     % store g2
    
        progressbar(jj/n_zone);
    end
end

%% VIS ----------------------------------------------------
idx_tau_vis=1;

% % g2(0) distribuion on halo
figure;
ax=tight_subplot(3,1);
for ii=1:3
    I=imagesc(ax(ii),squeeze(g2_mode{idx_tau_vis}(:,:,ii))');
    cbar=colorbar(ax(ii));
    colormap('inferno');
end
