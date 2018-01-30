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

param_id=param_array(:,1);
nparam=size(params,1);      % number of unique param-set

% group shot-ids by exp-param
id_in_param=cell(1,nparam);
for ii=1:nparam
    id_in_param{ii}=param_id(Ipar==ii);
end


% decompose param-set
dim_param=size(params,2);       % dimension of param-set

par_iter=cell(1,dim_param);     % iterated par-vec in each dim of par-space
for ii=1:dim_param
    par_iter{ii}=unique(params(:,ii));
end
n_par_iter=cellfun(@(p) numel(p),par_iter);     % num searched in each dim

% build index table to sort params
%   param-vect [p_i,q_j,...] is indexed by (i,j,...)
par_tab=NaN(size(params));
for ii=1:dim_param
    [~,par_tab(:,ii)]=ismember(params(:,ii),par_iter{ii});
end

