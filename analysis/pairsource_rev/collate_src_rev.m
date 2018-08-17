%% Collate pairsource characterisation data
% for 90 deg

%% load data
% path to pre-processed data
%   processed data from distinguishable source analysis (halo centered)
path_proc='C:\Users\HE BEC\bell_data\pairsource\90\rev';

% vars to load
vars2load={'g2_amp','g2_amp_sdev','n_mocc','n_mocc_err'};

% parse directory
dirlist=dir(path_proc);     
dirnames={dirlist.name};    
isfiles=cellfun(@(s) is_file(fullfile(path_proc,s)),dirnames);
flist=dirnames(isfiles);        % files in directory (only have .mat files in proc dir)
nfiles=numel(flist);

% load data
for ii=1:nfiles
    Sdata(ii)=load(fullfile(path_proc,flist{ii}),vars2load{:});
end

g2_amp=[Sdata.g2_amp];
g2_amp_sdev=[Sdata.g2_amp_sdev];
n_mocc=[Sdata.n_mocc];
n_mocc_err=[Sdata.n_mocc_err];
