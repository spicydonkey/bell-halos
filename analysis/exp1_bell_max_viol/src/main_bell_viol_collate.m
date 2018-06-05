% Statistical analysis on collated dataset
%
% DKS
% 20180605


%% CONFIG
% path to postprocessed data
path_data='C:\Users\David\Documents\bell\exp1_bell_inequality\postproc';


%% load data
%%% parse data directory
dir_data=dir(path_data);
fname_all={dir_data.name}';

% get mat data
b_file=cellfun(@(s) is_file(fullfile(path_data,s)),fname_all);
fname_mat=fname_all(b_file);

b_mat=cellfun(@(s) strcmp(s(end-2:end),'mat'),fname_mat);
fname_mat=fname_mat(b_mat);
ndata=numel(fname_mat);

% get file id
f_id=cellfun(@(s) str2num(strtok(extractAfter(s,'_'),'_')),fname_mat);
[f_id_sort,Is]=sort(f_id);
fname_mat_sort=fname_mat(Is);


%%% load
S_data=cell(ndata,1);

% loads everthing in the mat file (only ~40 MB total)
for ii=1:ndata
    S_data{ii}=load(fullfile(path_data,fname_mat_sort{ii}));
end


%% Preprocess
% get scattered atoms k-vectors from all and collate
nshot_run=cellfun(@(s) size(s.k_par{1},1),S_data);
cumsum_nshot_run=[0,cumsum(nshot_run)'];

k_col=cell(sum(nshot_run),2);       % prellocate; mf=-1 is irrelevant

for ii=1:ndata
    k_col(cumsum_nshot_run(ii)+1:cumsum_nshot_run(ii+1),:)=S_data{ii}.k_par{1}(:,1:2);
end

% some stats
Nsc_counts=shotSize(k_col);     % tot number of detected counts

figure;
hold on;
histogram(Nsc_counts(:,1));
histogram(Nsc_counts(:,2));
xlabel('Tot sc counts detected');
ylabel('num of shots');


%% Analysis
%% 1. g2 for all data


%% 2. Uncert from bootstrapping


%% 3. Correlation coefficient