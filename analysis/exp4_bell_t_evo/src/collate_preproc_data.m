%% collate preprocessed data
% DKS
% 2018-10-25

clear all;

%% configs
% path to pre-processed data
fpath='C:\Users\HE BEC\Dropbox\phd\data\bell_epr_2018\proc\exp4_tevo\preproc';

% vars to load
vars2load={'par_T','k_par'};

% parse directory
dirlist=dir(fpath);     
dirnames={dirlist.name};    
isfiles=cellfun(@(s) is_file(fullfile(fpath,s)),dirnames);
flist=dirnames(isfiles);        % files in directory (only have .mat files in proc dir)
nfiles=numel(flist);

%% load data
for ii=1:nfiles
    S(ii)=load(fullfile(fpath,flist{ii}),vars2load{:});
end

%% tidy data and naming
tau=cat(1,S.par_T);
k_tau=cat(1,S.k_par);

clearvars S;