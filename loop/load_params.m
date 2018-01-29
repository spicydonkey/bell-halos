%Analyse param log
%
%   params: unique param-sets scanned in experiment
%   id_in_param: shot-IDs categorised by param-set
%
% DK Shin

path_log = 'C:\Users\HE BEC\exp-data\bell\loop_tscan\LOG_parameters.txt';

param_log = load_logfile(path_log);
param_array = paramlog2array(param_log);

% get unique param-vecs and tag each shot with param-ID
[params,~,Ipar] = unique(param_array(:,2:end),'rows');   

id_shot=param_array(:,1);
nparam=size(params,1);      % number of unique param-set

% group shot-ids by exp-param
id_in_param=cell(1,nparam);
for ii=1:nparam
    id_in_param{ii}=id_shot(Ipar==ii);
end

% decompose param-set
mf_src=params(:,1);
t_raman=params(:,2);