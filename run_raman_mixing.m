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


%% Analyse count distribution
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

%% evaluate Rabi population
% NOTE: we separate the analysis here since mF=0 has a large background
cfg_bell_lo=[config_list.bell_lo];
sourcemf=[cfg_bell_lo.sourcemf];
ampRaman=[cfg_bell_lo.ampRaman];

P_rabi=cell(2,1);        % probability of population in mF=1 (UP) state
% cell1 - mF=0
% cell2 - mF=1

idx_mf=cellfun(@(x)find(sourcemf==x),{0,1},'UniformOutput',false);  % array index to mF=0/1
ampRaman_mf=cellfun(@(x)ampRaman(x),idx_mf,'UniformOutput',false);  % Raman amplitude separated for mF
for ii=1:2      % loop through mf=0/1
    P_rabi{ii}=cell(numel(idx_mf{ii}),1);   % preallocate cell-array
    
    for jj=1:numel(idx_mf{ii})
        idx_this=idx_mf{ii}(jj);
        
        % evaluate mean
        this_nn_halo_mean=cellfun(@(x)mean(x,3),nn_halo{idx_this},'UniformOutput',false);
        
        
        this_P_rabi=this_nn_halo_mean{2}./(this_nn_halo_mean{1}+this_nn_halo_mean{2});
        % 1 when all state in mf=1; 0 when all states in mf=0
        
        P_rabi{ii}{jj}=this_P_rabi;
    end
    
    P_rabi{ii}=cat(3,P_rabi{ii}{:});    % cell to matrix
end

%%% deal with large background for mf=0 case

%%% plot
%%% sort P_rabi(amp) according to ampRaman_mf
hfig_rabi_flop=figure(20);
for mm=1:2
    for ii=1:floor(14/2):14
        for jj=1:floor(19/2):19
            plot(ampRaman_mf{mm},squeeze(P_rabi{mm}(ii,jj,:)),'o',...
                'DisplayName',sprintf('%d: (%.2g,%.2g)',mm-1,ii,jj));
            hold on;
        end
    end
end
leg=legend('show','Location','eastoutside');
% set(leg,'Title','hello');
box on;
hold off;
xlabel('Raman Amplitude');
ylabel('$P(\uparrow)$');