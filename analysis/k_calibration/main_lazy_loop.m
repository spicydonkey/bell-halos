%% Lazy script for looping sensitivity analysis
% DKS
% 20181015

%%% WARNING: MUST NOT USE SAME VAR AS main_kcent_sensitivity.m which is
%%% iteratively called as a script
flag_load_override=true;
flag_save_data=true;


PATH_DIR='C:\Users\HE BEC\Documents\lab\bell_momentumspin\bell_epr_2018\proc\exp5_bell_yrot';
D=dir(PATH_DIR);
tfnames={D.name};
b_isfile=cellfun(@(s) is_file(fullfile(PATH_DIR,s)), tfnames);
filename=tfnames(b_isfile);
b_ismat=cellfun(@(s) strcmp(fileext(s),'.mat'),filename);
dataname=filename(b_ismat);
n_data=numel(dataname);

for II=2:n_data     % SKIP exp1_collated
% for II=2      % DEBUGGING
    %%%%%%%%%% LOAD SOME DATA
    if ~is_file(fullfile(PATH_DIR,dataname{II}))
        continue
    end
    tdataname=dataname{II};
    load(fullfile(PATH_DIR,tdataname));
    %%%%%%%%%% DATA IS LOADED FOR CALLING SCRIPT
    
    run('main_kcent_sensitivity.m');
end