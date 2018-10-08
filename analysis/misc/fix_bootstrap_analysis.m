%% Fix error in bootstrapping for SE of mean of correlator
%
% DKS
% 20181008

% config
path_dir='C:\Users\HE BEC\Documents\lab\bell_momentumspin\bell_epr_2018\proc\exp5_bell_yrot\bell_signal';
old_data='bell_20181008_oldstat.mat';
path_old_data=fullfile(path_dir,old_data);

% load data
load(path_old_data);    % load into workspace

% correct SE estimation by bootstrapping (scaled by sampling fraction)
for ii=1:numel(S)
    t_se_bs_old=S(ii).E_bootstrap_sdev;     % unscaled error estimate
    t_samp_frac=S(ii).n_frac_samp;          % bs sampling fraction
    
    S(ii).E_bootstrap_sdev_old=t_se_bs_old;         % store the old error
    S(ii).E_bootstrap_sdev=sqrt(t_samp_frac)*t_se_bs_old;    % sample size corrected
end