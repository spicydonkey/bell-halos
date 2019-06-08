%% SOM ANALYSIS: halo-integrated g2
% DKS
% 2019-03-20


%% CONFIGS setup
% data file
%   contains spin,momentum resolved vectors
fdata='C:\Users\HE BEC\Dropbox\PhD\projects\halo_metrology\analysis\exp4_tevo\exp4_20181029.mat';

% g2 
configs.g2.dk_lim=[-0.2,0.2];       % 0.2*[-1,1]
configs.g2.dk_n=15;                 % 15    

% vis -----------------------------------------------------------
config_fig = loadFigureConfig;      % load template


%% configure g2
dk_ed_vec=linspace(configs.g2.dk_lim(1),configs.g2.dk_lim(2),configs.g2.dk_n+1);
dk_cent_vec=edge2cent(dk_ed_vec);
[~,idx_dk0]=min(abs(dk_cent_vec));    % idx bin nearest zero
dk_ed={dk_ed_vec,dk_ed_vec,dk_ed_vec};
dk_grid_size=cellfun(@(k) length(k)-1,dk_ed);


%% load data
load(fdata);


%% debug: reduce data
% warning('Reducing data for speed.');
% for ii=1:numel(k_tau)
%     k_tau{ii}=k_tau{ii}(1:round(length(k_tau{ii})/10),:);
% end
 

%% MAIN ------------------------------------------------------
n_tau=length(k_tau);

g2_tau = cell(n_tau,1);
dk_tau = cell(n_tau,1);
g2mdl_tau = cell(n_tau,1);
h_tau = cell(n_tau,1);

for ii=1:n_tau
    t_K = k_tau{ii};
    [g2_tau{ii},dk_tau{ii},g2mdl_tau{ii},h_tau{ii}]=summary_disthalo_g2(t_K,dk_ed,0,1,1,1);
    % NOTE: manual flag for smoothing is INSIDE summary_disthalo_g2
end


%% VIS --------------------------------------------------------
% TODO
% comparison against smoothed(?) 
idx_tau_vis = 1;     % tau to disp


%% 3D vis: iso-surface
idx_spin_vis=2;     % display mJ=0,0 pair

t_g2_max = cellfun(@(x) max(x(:)),g2_tau{idx_tau_vis});
t_g2_iso_vis=linspace(1,t_g2_max(idx_spin_vis),5);
t_g2_iso_vis=t_g2_iso_vis(2:end-1);

figure;
hold on;
for gg=t_g2_iso_vis
    t_p = patch(isosurface(dk_tau{idx_tau_vis}{1},dk_tau{idx_tau_vis}{2},dk_tau{idx_tau_vis}{3},g2_tau{idx_tau_vis}{idx_spin_vis},gg),'FaceAlpha',0.1,'EdgeColor','none');
end
axis equal;
view(3);


%% 2D slice images: compare spins and tau
% config
idx_tau_vis = 6;     % tau to disp
ss_str = {'$\uparrow\uparrow$','$\downarrow\downarrow$','$\uparrow\downarrow$'};

figure('Name','g2_kint_2d','Units',config_fig.units,'Position',[0 0 12.8 3.75],'Renderer',config_fig.rend);
ax=tight_subplot(1,3,[0.1,0.1],[0.1,0.1],[0.1,0.1]);
for ii=1:3
    axes(ax(ii));
    I=imagesc('XData',dk_ed{1},'YData',dk_ed{2},'CData',squeeze(g2_tau{idx_tau_vis}{ii}(:,:,idx_dk0)));
    ax_pos=ax(ii).Position;
    if ii==3
        cbar=colorbar(ax(ii));
        title(cbar,'g^{(2)}')
    end
    ax(ii).CLim=[0,max(t_g2_max)];
    set(ax(ii),'Position',ax_pos);
    title(ss_str{ii});
    
    colormap('viridis');
    
    axis tight;
    axis equal;
    box on;
    set(ax(ii),'Layer','Top');
    ax(ii).FontSize=config_fig.ax_fontsize;
    ax(ii).LineWidth=config_fig.ax_lwid;
    ax(ii).XTickLabelMode='auto';
    ax(ii).YTickLabelMode='auto';
end

%% VIS: g2 surf
% data processing
g2_corr = cellfun(@(g2) 0.5*(g2{1} + g2{2}),g2_tau,'uni',0);     
g2_anti = cellfun(@(g2) g2{3},g2_tau,'uni',0);                      
g2_corrtyp = {g2_corr,g2_anti};

% config
ax_cmap={'red','blue'};
idx_tau_vis=1:6;

z_labels = {'$(\bar{g}^{(2)}_{\uparrow\uparrow} + \bar{g}^{(2)}_{\downarrow\downarrow})/2$','$\bar{g}^{(2)}_{\uparrow\downarrow}$'};

% xy-grid
[dk_2d_x,dk_2d_y]=ndgrid(dk_cent_vec,dk_cent_vec);

% vis
for ii=idx_tau_vis
    H=figure('Name','g2_kint','Units',config_fig.units,'Position',[0 0 4.3 8],'Renderer',config_fig.rend);
    H.Name=strcat(H.Name,'_',num2str(ii));
    
    ax=tight_subplot(2,1,[0.15,0.15],[0.1 0.1],[0.25 0.2]);
    for jj=1:2
        axes(ax(jj));
        
        ss = surf(gca,'XData',dk_2d_x,'YData',dk_2d_y,'ZData',squeeze(g2_corrtyp{jj}{ii}(:,:,idx_dk0)));
        
        colormap(gca,ax_cmap{jj});
        xlabel('$\Delta k_z$');
        ylabel('$\Delta k_x$');
%         zlabel('$\bar{g}^{(2)}$');
        zlabel(z_labels{jj});
        
        set(gca,'FontSize',config_fig.ax_fontsize);
        set(gca,'LineWidth',config_fig.ax_lwid);
        
        if jj==1
            titlestr = sprintf('%s = %0.3g ms','$\tau$',tau(ii));
            title(titlestr);
        end
    end
    
%     print_pnghr(H);       % save figs
end