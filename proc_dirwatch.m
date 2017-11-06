% groups exp-data (shots) according to same parameters assigned by
% dir_watch code
%
% PARAMS: PARAM1, PARAM2, ..., PARAMN
%
% NB to user: this script must be run while Matlab's pwd is configured to
% the directory containing data and the parameter log
%

logfile='LOG_parameters.txt';
paramlog=load_logfile(logfile);

dfiles=dir('d*.txt');
dfiles_id=cellfun(@(x) dfile2id(x),{dfiles(:).name}');  % get ordered shotID array

params=vertcat(paramlog(:).params);     % params as array
params=str2double(params);

params_list=unique(params,'rows');
[nparams_list,nparams]=size(params_list);

% categorise the dfiles into same-parameter bins
for ii=1:nparams_list
    tparam=params_list(ii,:);
    tdir=sprintf(repmat('%0.4g_',[1,nparams]),tparam);
    % directory name is "PARAM1_PARAM2_...PARAMN_"
    mkdir(tdir);
    
    % find files belonging to this param
    tbool=prod(params==tparam,2);
    tbool=(tbool==1);
    
    tid=vertcat(paramlog(tbool).id);
    
    % get those shots
    t_dbool=ismember(dfiles_id,tid);    
    t_dfiles={dfiles(t_dbool).name}';
    nfiles=length(t_dfiles);
    
    % move the files
    for jj=1:nfiles
        fprintf('%d/%d (%d/%d)\n',ii,nparams_list,jj,nfiles);

        % copyfile is very slow. FileRename is ~10x faster.
        %         copyfile(t_dfiles{jj},fullfile(tdir,t_dfiles{jj}));
        FileRename(t_dfiles{jj},fullfile(tdir,t_dfiles{jj}),'forced');
    end
end