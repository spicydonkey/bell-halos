%% Update sensitivity plot - g2 fitted amplitude
% DKS
% 20181018

%% configs
path_dir='C:\Users\HE BEC\Dropbox\phd\thesis\projects\bell_epr_steering\figS\k_calibration\g2_amp_fit\v1';

D=dir(path_dir);
tfnames={D.name};
b_isfile=cellfun(@(s) is_file(fullfile(path_dir,s)), tfnames);
filename=tfnames(b_isfile);
b_isfig=cellfun(@(s) strcmp(fileext(s),'.fig'),filename);
figname=filename(b_isfig);

nfigs=numel(figname);

%% figure makeup
% better y-label
H=cell(nfigs,1);
for ii=1:nfigs
    %%
    h=openfig(fullfile(path_dir,figname{ii}));
    figure(h);
    H{ii}=h;
    a=gca;
    
    ylabel('Norm. $g^{(2)}$ amplitude');
end

%% save data
if false
    %% Manually run this section
    path_save_dir='C:\Users\HE BEC\Dropbox\phd\thesis\projects\bell_epr_steering\figS\k_calibration\g2_amp_fit\v1_imp';
    for ii=1:nfigs
        fpath=fullfile(path_save_dir,figname{ii});
        
        saveas(H{ii},strcat(fpath,'.fig'),'fig');
        print(H{ii},strcat(fpath,'.svg'),'-dsvg');
    end
end