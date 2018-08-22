%% Create 2D images of Bell data (publication fig.1)
% DKS
% 20180822
%
% Prereq: prepare Nx3 "xyz"-vector array of data
%
flag_dbg=false;
%% load data
path_dir='C:\Users\HE BEC\Documents\lab\bell_momentumspin\bell_epr_2018\proc\fig1_halo_img';
path_data=fullfile(path_dir,'theta_pi2.mat');   % pi/2 rotation (ROT): too many shots
% path_data=fullfile(path_dir,'theta_0.mat');     % no rotation (SRC)

load(path_data);

% # counts in SRC:ROT = 825724 : 7142794
n_counts_src=825724;
% set counts same
xyz_orig=xyz;       % save orig
xyz=xyz(1:n_counts_src,:);      % get same num of counts

[n_v,v_ed,xyz_filt]=plot_pretty_halo(xyz,flag_dbg);