%% Lazy script for looping sensitivity analysis
% DKS
% 20181015

%%% WARNING: MUST NOT USE SAME VAR AS main_kcent_sensitivity.m which is
%%% iteratively called as a script
flag_load_override=true;
flag_save_data=true;


PATH_DIR='C:\Users\HE BEC\Documents\lab\bell_momentumspin\bell_epr_2018\proc\exp5_bell_yrot';
D=dir(PATH_DIR);

for II=1:numel(D)
    %%%%%%%%%% LOAD SOME DATA
    if ~is_file(fullfile(PATH_DIR,D(II).name))
        continue
    end
    load(fullfile(PATH_DIR,D(II).name));
    %%%%%%%%%% DATA IS LOADED FOR CALLING SCRIPT
    
    run('main_kcent_sensitivity.m');
end