% Charaterise single-beam Raman beam as a local operation
%
%

clear all; clc; close all;

%% configs
verbose=0;

% config files for characterising mixing
path_config='C:\Users\HE BEC\Documents\MATLAB\bell-halos\config\config_bell_mf_*';
config_files=dir(path_config);
config_files={config_files.name};   % name of config files


%% 
nconfigs=numel(config_files);
config_list=[];

nn_halo=cell(nconfigs,1);
azim=cell(nconfigs,1);
elev=cell(nconfigs,1);

clearvars config_list;

for ii=1:nconfigs
    % load the config "configs"
    run(config_files{ii});
    config_list(ii)=configs;    % store config
    
    % run analysis
    [nn_halo{ii},azim{ii},elev{ii}]=sph_zone_analysis(configs,verbose);
end