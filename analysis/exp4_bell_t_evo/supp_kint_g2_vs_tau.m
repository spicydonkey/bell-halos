%% SOM ANALYSIS: halo-integrated g2
% DKS
% 2019-03-20


%% CONFIGS setup
% data file
%   contains spin,momentum resolved vectors
fdata='C:\Users\HE BEC\Dropbox\PhD\projects\halo_metrology\analysis\exp4_tevo\exp4_20181029.mat';

% g2 
configs.g2.dk_lim=[-0.1,0.1];
configs.g2.dk_n=9;

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
end

%% VIS --------------------------------------------------------
% TODO