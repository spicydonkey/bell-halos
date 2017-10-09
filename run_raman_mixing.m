% Charaterise single-beam Raman beam as a local operation
%
%

clear all; clc; close all;

%% CONFIG
override_config=1;
    load_config_default=0;
        config_default_id=1;
    VERBOSE=0;

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
    
    verbose=configs.flags.verbose;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OVERRIDE CONFIG DEFINED FROM FILE
    if override_config
        verbose=VERBOSE;
        
        configs.zone.nazim=100;
        configs.zone.nelev=50;
        
        if load_config_default
            run(sprintf('config_default_%d',config_default_id));
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % store config
    config_list(ii)=configs;    
    
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
% TODO

%%% sort Rabi cycle data wrt Raman amplitude
for ii=1:2
    % ascending sort Raman amplitude vector; get sorting index
    [ampRaman_mf{ii},this_Isort]=sort(ampRaman_mf{ii},'ascend');
    
    % sort the Rabi cycle
    P_rabi{ii}=P_rabi{ii}(:,:,this_Isort);
end


%% Evaluate rotation angle from the two-level population ratio
% formula: Pr(ORIGIGNAL_STATE)=cos(ROT_ANGLE/2)^2
%   ROT_ANGLE >= 0 by convention
% TODO - can only tell rotation angle modulo pi

Theta=cellfun(@(Pr) 2*acos(sqrt(Pr)),P_rabi,'UniformOutput',false);
Theta{1}=pi-Theta{1};   % mf=0
% TODO - why not Theta{1}=pi+Theta{1} ??

%% plot Rabi cycle
hfig_rabi_flop=figure(20);

% same spherical grid for each experiment!
nazim=size(azim{1},2)-1;    % subtract 1: bin centers from edge
nelev=size(azim{1},1)-1;

% grid locations to show Rabi flopping
nazim_div=5;
nelev_div=5;
% nazim_div=floor(nazim/1);
% nelev_div=floor(nelev/1);

% misc
pcolors=distinguishable_colors(nazim_div*nazim_div);
pmarkers={'o','^'};
plinestyles={'-','--'};

for mm=1:2
    counter=1;
    for ii=1:ceil(nelev/nelev_div):nelev
        for jj=1:ceil(nazim/nazim_div):nazim
            plot(ampRaman_mf{mm},squeeze(P_rabi{mm}(ii,jj,:)),...
                'LineStyle',plinestyles{mm},'Color',pcolors(counter,:),...
                'Marker',pmarkers{mm},...   %'MarkerFaceColor',pcolors(counter,:),
                'DisplayName',sprintf('%d: (%.2g,%.2g)',mm-1,ii,jj));
            hold on;
            counter=counter+1;
        end
    end
end
leg=legend('show','Location','eastoutside');
title(leg,'$m_F$: (elev ii, azim jj)');
box on;
hold off;
xlabel('Raman Amplitude');
ylabel('$P(\uparrow)$');

%% plot rotation angle
hfig_theta_lo=figure(21);

% grid locations to show Rabi flopping
nazim_div=5;
nelev_div=5;
% nazim_div=floor(nazim/1);
% nelev_div=floor(nelev/1);

for mm=1:2
    counter=1;
    for ii=1:ceil(nelev/nelev_div):nelev
        for jj=1:ceil(nazim/nazim_div):nazim
            plot(ampRaman_mf{mm},squeeze(Theta{mm}(ii,jj,:)),...
                'LineStyle',plinestyles{mm},'Color',pcolors(counter,:),...
                'Marker',pmarkers{mm},...   %'MarkerFaceColor',pcolors(counter,:),
                'DisplayName',sprintf('%d: (%.2g,%.2g)',mm-1,ii,jj));
            hold on;
            counter=counter+1;
        end
    end
end
leg=legend('show','Location','eastoutside');
title(leg,'$m_F$: (elev ii, azim jj)');
box on;
hold off;
xlabel('Raman Amplitude');
ylabel('$\Theta$');

%% Plot rotation map
lon=rad2deg(azim{1});
lat=rad2deg(elev{2});

numElevEdgesNaN=11;

hFiltGauss=fspecial('gaussian');

%%% ECKERT4 projection
for ii=1:size(Theta{2},3)
    figure;
    
    Th=Theta{2}(:,:,ii);
    
    % remove data near caps
    Th([1:numElevEdgesNaN,end-numElevEdgesNaN:end],:)=NaN;
    
    % apply gaussian filter
    Th=imfilter(Th,hFiltGauss);
    
    plotFlatMap(lat,lon,Th,'eckert4');

    % misc
    ax=gca;
    
    tstr=sprintf('$m_F=+1$, Raman Amp = %0.2g',ampRaman_mf{2}(ii));
    title(tstr);
    colormap('magma');
    hcb=colorbar('Southoutside');
    hcb.Label.String='$\psi(\theta,\phi)$';
    hcb.TickLabelInterpreter='latex';
    hcb.Label.Interpreter='latex';
    ax.FontSize=12;
    
    drawnow;
    
    % save
    fnamethis=sprintf('theta_eck_mf1_%d.png',ii);
    saveas(gcf,fnamethis);
end

%%% rectangle
for ii=1:size(Theta{2},3)
    figure;
    
    Th=Theta{2}(:,:,ii);
    
    % remove data near caps
    Th([1:numElevEdgesNaN,end-numElevEdgesNaN:end],:)=NaN;
    
    % apply gaussian filter
    Th=imfilter(Th,hFiltGauss);
    
    plotFlatMap(lat,lon,Th,'rect');    

    % misc
    ax=gca;
    
    tstr=sprintf('$m_F=+1$, Raman Amp = %0.2g',ampRaman_mf{2}(ii));
    title(tstr);
        
    xlim([-180,180]);
    ylim([-90,90]);
    
    colormap('magma');
    hcb=colorbar('Southoutside');
    hcb.Label.String='$\psi(\theta,\phi)$';
    hcb.TickLabelInterpreter='latex';
    hcb.Label.Interpreter='latex';
    ax.FontSize=12;
    
    drawnow;
    
    % save
    fnamethis=sprintf('theta_rect_mf1_%d.png',ii);
    saveas(gcf,fnamethis);
end