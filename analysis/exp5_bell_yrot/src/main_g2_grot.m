%% g2 analysed data from global rotation
%
% DKS
% 2018-08-14

%% load data
% path to pre-processed data
path_proc='C:\Users\HE BEC\Documents\lab\bell_momentumspin\bell_epr_2018\proc\exp5_bell_yrot';

% vars to load
vars2load={'par_T','g2mdl','g2','dk',...
    'g2_sdev',...
    'E_bootstrap_sdev',...
    'g2anti_par','g2corr_par','E_par'};

% parse directory
dirlist=dir(path_proc);     
dirnames={dirlist.name};    
isfiles=cellfun(@(s) is_file(fullfile(path_proc,s)),dirnames);
flist=dirnames(isfiles);        % files in directory (only have .mat files in proc dir)
nfiles=numel(flist);

% load data
% Sdata=cell(nfiles,1);   % preallocate
for ii=1:nfiles
    S(ii)=load(fullfile(path_proc,flist{ii}),vars2load{:});
end



