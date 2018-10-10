%% Create 2D images of halos mixing mJ
% DKS
% 20181010

%% load data
% need ZXY processed data
path_dir='C:\Users\HE BEC\bell_data';
fname_data='expS1_load.mat';
path_data=fullfile(path_dir,fname_data);

load(path_data);        % load data


%% Process data to vis
% distinguish ZXY data by scanned param (pulse tau)
%   xyz0 is x,y,z' where the z' is shifted to window center
xyz0_par = cell(nparam,1);
% z_shift=-4.0809*[0,0,mean(configs.load.window{1})];
z_shift=[0,0,-1.66];
for ii=1:nparam
    xyz0_par{ii}=circshift(cat(1,zxy{b_paramset(:,ii)}),-1,2)+z_shift;
end

%%% filter
% detector hot-spots - cylinrical filter
det_cent=[0,3e-3,0];
hotspot_rad=30e-3;
xyz0_filt_1=cellfun(@(x) cylindercull(x,det_cent,[hotspot_rad,Inf],3),xyz0_par,'UniformOutput',false);

% box window
crop_window={[-30e-3,30e-3],[-1.2e-2,1.5e-2],[-0.1,0.12]};
xyz0_filt_2 = cellfun(@(x) boxcull(x,crop_window),xyz0_filt_1,'UniformOutput',false);

%%% integrate along Y-dir
xz=cellfun(@(xyz) xyz(:,[1,3]),xyz0_filt_2,'UniformOutput',false);

%% number histogram
len_bin=0.05e-3;      % bin size
ed_bin=cellfun(@(x) linspace(x(1),x(2),round(diff(x)/len_bin)),crop_window,'UniformOutput',false);
ed_bin=ed_bin([1,3]);

N=cellfun(@(x) nhist(x,ed_bin),xz,'UniformOutput',false);   % tot atom # detected in bin

% avg atom number density
det_qe=0.1;
n_shot_par=sum(b_paramset,1);
vol_bin=1e9*diff(crop_window{2})*len_bin^2;     % vol of bin in mm^3

n=cell(nparam,1);
for ii=1:nparam
    n{ii}=N{ii}*(1/det_qe)*(1/n_shot_par(ii))*(1/vol_bin);
end

% Gaussian filter
gaussfilt_sig=10;
n_filt=cellfun(@(x) imgaussfilt(x,gaussfilt_sig),n,'UniformOutput',false);


%% vis
hfig=cell(nparam,1);
for ii=1:nparam
    %%
    hfig{ii}=figure('Units', 'normalized', 'Position', [0.2,0.2,0.2,0.6]);
    hfig{ii}.Renderer='painters';
    
    bgAxes=axes('Position',[0,0,1,1],'XColor','none','YColor','none',...
        'XLim',[0,1],'YLim',[0,1]);     % bg axes
    
    text_mj={'1','0','-1'};
    pos_x=0.10;
    pos_y=[0.24,0.48,0.72];
    for jj=1:3
        text(pos_x,pos_y(jj),text_mj{jj},'FontSize',16,'Rotation',90,...
            'HorizontalAlignment','center','VerticalAlignment','middle',...
            'FontWeight','bold');
    end
    text(pos_x-0.075,pos_y(2),'$m_J$','FontSize',16,'Rotation',90,...
        'HorizontalAlignment','center','VerticalAlignment','middle');
    
    ax_img=axes('Position',[0.1,0.1,1,0.8],'XColor','none','YColor','none',...
        'XLim',[0,1],'YLim',[0,1]);
    pimg=imagesc(1e3*ed_bin{1},1e3*ed_bin{2},log10(n_filt{ii})');    

    
    get(ax_img);
    ax_img.YDir='normal';
    axis equal;
    axis tight;
    ax_img.FontSize=12;
    xlabel('$x$ [mm]');
    ylabel('$z$ [mm]');
    title(strcat(sprintf('%0.2g',1e6*par_T(ii)),' $\mu$s'));
    colormap('inferno');
    cbar=colorbar('FontSize',11);
    cbar.Label.Interpreter='latex';
    cbar.Label.String='$\log($ avg atom number density $[\textrm{mm}^-3]$ $)$';
    caxis_lim=[-2,1];
    caxis(caxis_lim);
    cbar.Ticks=caxis_lim(1):1:caxis_lim(2);
end

%% save data
if false        % must be called manually
    %% Manually run this section
    path_save_dir='C:\Users\HE BEC\Dropbox\PhD\thesis\projects\bell_epr_steering\figS\img_halo_mix';
    
    for ii=1:nparam
        fname=sprintf('fig_%d_%s',ii,getdatetimestr);
        fpath=fullfile(path_save_dir,fname);
        
%         saveas(hfig{ii},fpath,'fig');     % HUGE!
        print(hfig{ii},fpath,'-dsvg');
    end
end