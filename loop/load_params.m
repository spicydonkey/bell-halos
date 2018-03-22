%Analyse param log
%
%   params: unique param-sets scanned in experiment
%   id_in_param: shot-IDs categorised by param-set
%
%   TODO:
%       [ ] last step is specific to param scan -- need to generalise
%
% DK Shin
%

%% CONFIG
% dir_base='X:\expdata\spinmom_bell\loop_tscan';

dir_base='X:\expdata\spinmom_bell\loop_tscan3';
path_log=fullfile(dir_base,'raw','temp_LOG_parameters.txt');


%% main
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
% deal with single-param scan
if numel(n_par_iter)==1
    n_par_iter=[n_par_iter,1];
end

% build index table to sort params
%   param-vect [p_i,q_j,...] is indexed by (i,j,...)
par_tab=NaN(size(params));
for ii=1:dim_param
    [~,par_tab(:,ii)]=ismember(params(:,ii),par_iter{ii});
end


% build reverse-table
%   data location at (i,j,...) retrieves param-vect
parvec_tab=cell(n_par_iter);
for ii=1:nparam
    this_idx=num2cell(par_tab(ii,:));     % location in table 
    parvec_tab{this_idx{:}}=params(ii,:);
end


%% Scan dependent
% summary of experimental params
% 1. 201801XX: scan 1,2
% mf=par_iter{1};
% traman=par_iter{2};

% 2. 20180222: scan 3
traman=par_iter{1};