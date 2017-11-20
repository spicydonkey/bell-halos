% Load Raman data

% User config
%%% Bell data
% path_config='config_bell_load_3ms_0K38_medn_1.m';
% path_config='config_bell_load_3ms_0K38_medn_2.m';
% path_config='config_bell_load_3ms_0K29_1.m';
% path_config='config_bell_load_3ms_0K29_2.m';
% TODO
path_config='config_bell_load_3ms_0K38_lown_1.m';

%%% loop data
% path_config='config_loop_load_3ms_0K382.m';
% path_config='config_loop_load_3ms_0K2889.m';

% vars to save to output
vars_save={'configs','path_config',...
    'fullpath_config',...
    'fname_save','path_save',...
    'bec_cent',...
    'Nsc_mean','Nsc_std',...
    'dk_mean','dk_std',...
    'halo_k',...
    };

%% Main
t_main_start=tic;

% configure
dirsource=fileparts(mfilename('fullpath'));     % device independent source directory
fullpath_config=fullfile(dirsource,'config/load',path_config);   % full path to configs dir
run(fullpath_config);

% update configuration
configs=updateConfig(configs);

% set up this run's ID and misc paths
run_id=getdatetimestr;
fname_save=[mfilename,'__',run_id];
path_save=fullfile(fileparts(configs.files.path),'arch',fname_save);
if configs.flags.savedata || configs.flags.savefigs
    mkdir(path_save);   % create dir to output any files to save to disk
end

%% get ZXY counts
dpath=fullfile(configs.files.path,'d');
[txy,fout]=load_txy(dpath,configs.load.id,configs.load.window,...
    configs.load.mincount,configs.load.maxcount,...
    configs.load.rot_angle,1,...
    configs.flags.verbose,configs.flags.graphics);
zxy=txy2zxy(txy);
nshots=size(zxy,1);

%% get halos for mf=0,1
% preallocate vars
halo_k=cell(nshots,2);
Nsc_mean=NaN(1,2);
Nsc_std=NaN(1,2);
dk_mean=NaN(1,2);
dk_std=NaN(1,2);
bec_cent=cell(1,2);
hfig_halos=cell(1,2);

for jj=1:2
    [halo_k(:,jj),Nsc_mean(jj),Nsc_std(jj),dk_mean(jj),dk_std(jj),bec_cent{jj},hfig_halos{jj}]=halo_2bec(zxy,...
        configs.bec{jj}.pos{1},configs.bec{jj}.pos{2},...
        configs.bec{jj}.Rmax,configs.bec{jj}.r_th,...
        configs.halo{jj}.dR,configs.halo{jj}.elev_max,...
        configs.flags.verbose);
end

clearvars txy zxy;   % clean workspace


%% Save results
if configs.flags.savedata
    % TODO - package this into a function
    % TODO - if datafile already exists, MOVE it
    
    % parse list of vars user wants to save and check against workspace if exists
    varsExist=cell(size(vars_save));
    varCounter=0;
    for i = 1:length(vars_save)
        tVar=vars_save{i};
        if ~exist(tVar,'var')
            warning(['Variable "',tVar,'" does not exist.']);
            continue;
        end
        varCounter=varCounter+1;
        varsExist{varCounter}=tVar;
    end
    varsExist=varsExist(1:varCounter);
    
    dir_save=fileparts(path_save);
    if ~exist(dir_save,'dir')
        warning('Directory to save data %s does not exist. Creating directory.',dir_save);
        mkdir(dir_save);
    end
    save([path_save,'.mat'],varsExist{:});
end

%% END
t_main_end=toc(t_main_start);
disp('-----------------------------------------------');
fprintf('Total elapsed time (s): %7.1f\n',t_main_end);
disp('===================ALL TASKS COMPLETED===================');
