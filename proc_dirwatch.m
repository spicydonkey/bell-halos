% organise data into shots with same experimental parameter
%
% PARAMS: RAMAN_AMPLITUDE
%
% TO USER: go to the directory containing the data
%

logfile='LOG_parameters.txt';
paramlog=load_logfile(logfile);

dfiles=dir('d*.txt');
dfiles_id=cellfun(@(x) dfile2id(x),{dfiles(:).name}');

params=vertcat(paramlog(:).params);
params=str2double(params);

params_list=unique(params,'rows');
[nparams_list,nparams]=size(params_list);

for ii=1:nparams_list
    tparam=params_list(ii,:);
    tdir=sprintf(repmat('%0.4g_',[1,nparams]),tparam);
    mkdir(tdir);
    
    tbool=prod(params==tparam,2);
    tbool=(tbool==1);
    
    tid=vertcat(paramlog(tbool).id);
    
    t_dbool=ismember(dfiles_id,tid);
    t_dfiles={dfiles(t_dbool).name}';
    nfiles=length(t_dfiles);
    
    for jj=1:nfiles
        fprintf('%d/%d (%d/%d)\n',ii,nparams_list,jj,nfiles);
%         copyfile(t_dfiles{jj},fullfile(tdir,t_dfiles{jj}));
          FileRename(t_dfiles{jj},fullfile(tdir,t_dfiles{jj}),'forced');  
    end
end