% Bell test numerical simulation

% clear all; close all; clc;

OVERRIDE_CONFIG_FLAG=true;

%% general config
% datetimestr=datestr(datetime,'yyyymmdd_HHMMSS');    % timestamp when function called
datetimestr='test';
path_data_dir='C:\Users\HE BEC\bell\temp';
path_data_bell=sprintf('bell_sim_%s',datetimestr);
path_data_loop=sprintf('loop_sim_%s',datetimestr);


%% BELL TEST
%%% experiment
n_shot=1e3;             % number of shots
n_corr_pairs=50*ones(n_shot,1);        % number of correlated pairs generated from source
det_qe=1*ones(1,2);       % detection efficiency

% momentum width (implemented on scattered)
dk_dither_sd=zeros(1,3);
% dk_dither_sd=3e-2*[1,1,1];

% local operation
% TODO - we can even load the phase map from experiment!
% fun_localoper=@(th,phi) pi*sin(sqrt(76)*th/(2*pi)).*sin(sqrt(7)*phi*(2*pi)/(pi/2));
% fun_localoper=@(th,phi) pi*sin(20*th).*sin(20*phi);
% fun_localoper=@(th,phi) 1.5*phi+pi/2;
fun_localoper=@(th,phi) 2*phi;
% fun_localoper=@(th,phi) 0*phi;              % no local operation

%%% HALO
configs.halo{1}.string='$m_F=0$';
configs.halo{2}.string='$m_F=1$';

% boosts in this case misalign momentum correlator
% configs.halo{1}.boost=-1e-2/sqrt(3)*[1,1,1];
% configs.halo{2}.boost=1e-2/sqrt(3)*[1,1,1];
configs.halo{1}.boost=[0,0,0];
configs.halo{2}.boost=[0,0,0];

%%% Spherical zones
configs.zone.nazim=100;
configs.zone.nelev=50;

configs.zone.binmethod=1;
configs.zone.binwidth=0.05;

run('gen_ideal_corr_source.m');

% save outputs
% required variables for analysis currently are:
%   azim_grid elev_grid E 
varstosave={'azim_grid','elev_grid','E'};
save(fullfile(path_data_dir,path_data_bell),varstosave{:});

%% MIXING PULSE
% configure:
%   
verbose=1;
nazim=100;
nelev=50;

run('gen_ideal_mixing.m');

% save outputs
% required variables for analysis currently are:
%   ampRaman_mf azim elev Theta
varstosave={'ampRaman_mf','azim','elev','Theta'};
save(fullfile(path_data_dir,path_data_loop),varstosave{:});


%% ANALYSIS
raman_amp=NaN;      % for simulated local mixing

run('run_bell.m');



%% clear override
clearvars OVERRIDE_CONFIG_FLAG;     % REMOVE OVERRIDE FLAG WHEN DONE