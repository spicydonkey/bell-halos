%% Collate all runs at maximal Bell signal
%
% DKS
% 20180815


%% load data
% path to pre-processed data
path_proc='C:\Users\HE BEC\Documents\lab\bell_momentumspin\bell_epr_2018\proc\exp1_bell_max_viol';

% vars to load
vars2load={'par_T','k_par','nparam','n_mf','configs'};
%     'g2','dk','g2mdl','idx_dk0',...
%     'g2anti_par','g2corr_par','E_par',...
%     'g2_sdev','E_bootstrap_sdev'};

% parse directory
dirlist=dir(path_proc);     
dirnames={dirlist.name};    
isfiles=cellfun(@(s) is_file(fullfile(path_proc,s)),dirnames);
flist=dirnames(isfiles);        % files in directory (only have .mat files in proc dir)
nfiles=numel(flist);

% load data
% Sdata=cell(nfiles,1);   % preallocate
for ii=1:nfiles
    S_collated(ii)=load(fullfile(path_proc,flist{ii}),vars2load{:});
end


%% Collate data
% k_par
k_par={};
for ii=1:nfiles
    k_par=vertcat(k_par,S_collated(ii).k_par{1});
end
k_par={k_par};      % format k_par as a 1-by-nparam cell-array

% par_T
par_T=S_collated(1).par_T;   % all par_T in each run should be identical (5us)

% nparam
nparam=S_collated(1).nparam; % same as above and should be 1

% n_mf
n_mf=S_collated(1).n_mf;     % same as above
