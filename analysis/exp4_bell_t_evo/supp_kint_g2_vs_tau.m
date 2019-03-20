%% SOM ANALYSIS: halo-integrated g2
% DKS
% 2019-03-20


%% CONFIGS setup
% data file
%   contains spin,momentum resolved vectors
fdata='C:\Users\HE BEC\Dropbox\PhD\projects\halo_metrology\analysis\exp4_tevo\exp4_20181029.mat';

% g2 
configs.g2.dk_lim=[-0.2,0.2];       % 0.2*[-1,1]
configs.g2.dk_n=15;                 % 15    

% vis -----------------------------------------------------------
config_fig = loadFigureConfig;      % load template


%% configure g2
dk_ed_vec=linspace(configs.g2.dk_lim(1),configs.g2.dk_lim(2),configs.g2.dk_n+1);
dk_cent_vec=dk_ed_vec(1:end-1)+0.5*diff(dk_ed_vec);
[~,idx_dk0]=min(abs(dk_cent_vec));    % idx bin nearest zero
dk_ed={dk_ed_vec,dk_ed_vec,dk_ed_vec};
dk_grid_size=cellfun(@(k) length(k)-1,dk_ed);


%% load data
load(fdata);


%% debug: reduce data
% warning('Reducing data for speed.');
% for ii=1:numel(k_tau)
%     k_tau{ii}=k_tau{ii}(1:round(length(k_tau{ii})/10),:);
% end
 

%% MAIN ------------------------------------------------------
n_tau=length(k_tau);

g2_tau = cell(n_tau,1);
dk_tau = cell(n_tau,1);
g2mdl_tau = cell(n_tau,1);
h_tau = cell(n_tau,1);

for ii=1:n_tau
    t_K = k_tau{ii};
    [g2_tau{ii},dk_tau{ii},g2mdl_tau{ii},h_tau{ii}]=summary_disthalo_g2(t_K,dk_ed,0,1,1,1);
    % NOTE: manual flag for smoothing is INSIDE summary_disthalo_g2
end


%% VIS --------------------------------------------------------
% TODO
% comparison against smoothed(?) 
idx_tau_vis = 1;     % tau to disp


%% 3D vis: iso-surface
idx_spin_vis=2;     % display mJ=0,0 pair

figure;
hold on;
for gg=10:10:40
    t_p = patch(isosurface(dk_tau{idx_tau_vis}{1},dk_tau{idx_tau_vis}{2},dk_tau{idx_tau_vis}{3},g2_tau{idx_tau_vis}{idx_spin_vis},gg),'FaceAlpha',0.1,'EdgeColor','none');
end
axis equal;
view(3);


%% 2D slice images: compare spins and tau
figure; 
ax=tight_subplot(1,3);
for ii=1:3
    imagesc(ax(ii),squeeze(g2_tau{idx_tau_vis}{ii}(:,:,5)));
    t_cbar=colorbar(ax(ii));
    ax(ii).CLim=[0,60];
end
