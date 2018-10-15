%% edit sensitivity plot
% DKS
% 20181016

%% configs
path_dir='C:\Users\HE BEC\Dropbox\PhD\thesis\projects\bell_epr_steering\figS\k_calibration\v1';

D=dir(path_dir);
tfnames={D.name};
b_isfile=cellfun(@(s) is_file(fullfile(path_dir,s)), tfnames);
filename=tfnames(b_isfile);
b_isfig=cellfun(@(s) strcmp(fileext(s),'.fig'),filename);
figname=filename(b_isfig);

nfigs=numel(figname);

%% figure makeup
H=cell(nfigs,1);
TAU=NaN(nfigs,1);
for ii=1:nfigs
    %%
    % [ ] resize figure
    % [ ] annotate with tau
    h=openfig(fullfile(path_dir,figname{ii}));
    figure(h);
    H{ii}=h;
    a=gca;
    
    tau=str2double(h.Name(12:end));
    str_tau=sprintf('%s%0.3g %s','$\tau=$',tau,'$\mu$s');
    TAU(ii)=tau;
    
%     f.Position=[0.2 0.2 0.2 0.3];
%     a.Position=[0.15 0.15 1-0.3 1-0.3];
    
    text('Units','normalized','Position',[0.95 0.9],'String',str_tau,...
        'FontSize',14,'VerticalAlignment','middle','HorizontalAlignment','right');    
end

%% save data
if false
    %% Manually run this section
    path_save_dir='C:\Users\HE BEC\Dropbox\PhD\thesis\projects\bell_epr_steering\figS\k_calibration\v1_improved';
    for ii=1:nfigs
        fname=sprintf('fig_%0.2g_%s',TAU(ii),getdatetimestr);
        fpath=fullfile(path_save_dir,fname);
        
        saveas(H{ii},strcat(fpath,'.fig'),'fig');
        print(H{ii},strcat(fpath,'.svg'),'-dsvg');
    end
end