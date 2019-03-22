%% SOM ANALYSIS: localised g2
% DKS
% 2019-03-20

%% CONFIGS -------------------------------------------------
% data file
%   contains spin,momentum resolved vectors
fdata='C:\Users\HE BEC\Dropbox\PhD\projects\halo_metrology\analysis\exp4_tevo\exp4_20181029.mat';

% g2 
configs.g2.dk_lim=0.05*[-1,1];       
configs.g2.dk_n=1;                 

% k-modes
configs.mode.alpha=pi/10;           % cone half-angle
configs.mode.lim_az=pi*[0,1];       % limit of azim angle (exc +pi)
configs.mode.lim_el=pi/2*[-1,1];
configs.mode.n_az=30;               % equispaced bins
configs.mode.n_el=31;               % should be ODD

% post-processing
configs.post.el_trunc=deg2rad(50);           % elevation angle threshold to truncate (rad)

% vis 
config_fig = loadFigureConfig;      % load template

%% load data -------------------------------------------------
load(fdata);

%%% debug: reduce data
% warning('Reducing data for speed.');
% k_tau{1}=k_tau{1}(1:3000,:);
% for ii=1:numel(k_tau)
%     k_tau{ii}=k_tau{ii}(1:round(length(k_tau{ii})/10),:);
% end

% cull tau
% warning('culling Tau');
% k_tau=k_tau(1);

%% MAIN ------------------------------------------------------
% % configure g2 analysis
dk_ed_vec=linspace(configs.g2.dk_lim(1),configs.g2.dk_lim(2),configs.g2.dk_n+1);
dk_cent_vec=dk_ed_vec(1:end-1)+0.5*diff(dk_ed_vec);
[~,idx_dk0]=min(abs(dk_cent_vec));    % idx bin nearest zero
dk_ed={dk_ed_vec,dk_ed_vec,dk_ed_vec};
dk_grid_size=cellfun(@(k) length(k)-1,dk_ed);

% % configure k-modes
vaz=linspace(configs.mode.lim_az(1),configs.mode.lim_az(2),configs.mode.n_az+1);
vaz=vaz(1:end-1);       % exclude last
vel=linspace(configs.mode.lim_el(1),configs.mode.lim_el(2),configs.mode.n_el);
[~,idx_el0]=min(abs(vel));      % index to elev-angle elem closest to equator

[gaz,gel]=ndgrid(vaz,vel);    % AZ-EL grid
n_zone=numel(gaz);

% % localised g2 analysis
n_tau=length(k_tau);

g2_mode = cell(n_tau,1);

for ii=1:n_tau
    % preallocate
    g2_mode{ii}=NaN(configs.mode.n_az,configs.mode.n_el,3);

    t_k = k_tau{ii};    

    % g2 spatial distribution
    progressbar(0)
    for jj=1:n_zone
        [iaz,iel]=ind2sub(size(gaz),jj);
        taz=gaz(jj);
        tel=gel(jj);

        tk_mode = cellfun(@(x) inDoubleCone(x,taz,tel,configs.mode.alpha),t_k,'uni',0);

        % g2
        t_g2 = summary_disthalo_g2(tk_mode,dk_ed,0,0,0,0);

        g2_mode{ii}(iaz,iel,:) = [t_g2{:}];     % store g2
    
        progressbar(jj/n_zone);
    end
end


%% post processing
% s% truncate elevation angle limits
b_trunc=(abs(gel)>configs.post.el_trunc);    % boolean indicator to truncate

for ii=1:n_tau
    g2_mode{ii}(repmat(b_trunc,[1,1,3]))=NaN;
end


% % fill all (th,phi) by inversion symmetry
[gaz_filled,gel_filled,~]=autofill_cent_symm(gaz,gel,0*gaz);
g2_mode_filled=cell(size(g2_mode));
for ii=1:n_tau
    for jj=1:3 
        [~,~,g2_mode_filled{ii}(:,:,jj)] = autofill_cent_symm(gaz,gel,g2_mode{ii}(:,:,jj));
    end
end

vaz_filled = gaz_filled(:,1);       % azim angles 


% % data processing
g2_anti = cellfun(@(x) mean(x(:,:,1:2),3),g2_mode_filled,'uni',0);
g2_corr = cellfun(@(x) x(:,:,3),g2_mode_filled,'uni',0);


%% VERBOSE SUMMARY ----------------------------------------
g2_mode_avg = cellfun(@(x) arrayfun(@(ii) meanall(x(:,:,ii),'omitnan'),1:3),g2_mode,'uni',0);

fprintf('-------------------------------------\n');
for ii=1:n_tau
    fprintf('tau=%0.3g:\t%0.3g\t%0.3g\t%0.3g\n',tau(ii),g2_mode_avg{ii})
end
fprintf('-------------------------------------\n');


%% VIS ----------------------------------------------------
idx_tau_vis_dbg=1;      % for debugging test
idx_tau_vis=1:6;


%% VIS: localised g2 distribution
ss_str = {'$\uparrow\uparrow$','$\downarrow\downarrow$','$\uparrow\downarrow$'};

figure('Name','g2_kloc_all','Units',config_fig.units,'Position',[0 0 7 8.5],'Renderer',config_fig.rend);
ax=tight_subplot(3,1,[0,0],[0.12,0.10],[0.10,0.10]);
for ii=1:3
    axes(ax(ii));
    plotFlatMapWrappedRad(gaz_filled,gel_filled,squeeze(g2_mode_filled{idx_tau_vis_dbg}(:,:,ii)),'rect','texturemap');
    
    ax_pos0 = plotboxpos(ax(ii));
    cbar=colorbar();
    cbar.Position(3)=0.03;
    cbar.Position(1)=sum(ax_pos0([1,3]))+cbar.Position(3);
    cbar.Position([2,4])=ax_pos0([2,4]);
    
    if ii==1
        title(cbar,'g^{(2)}');
    end
    colormap('magma');
    set(gca,'Layer','Top');
    axis tight;
    axis equal;
    box on;
    ax(ii).FontSize=config_fig.ax_fontsize;
    ax(ii).LineWidth=config_fig.ax_lwid;
    
    if ii==3
        xlabel('$\theta$ (deg)');
    else
        set(gca,'XTickLabel','');
    end
    if ii==2
        ylabel('$\phi$ (deg)');
    end

    tt=text(0,90,ss_str{ii},'HorizontalAlignment','center','VerticalAlignment','top','BackgroundColor','w');
    uistack(tt,'bottom');
    
    hatch_axis(gca);
end

%% VIS: compare ANTI- vs CORR
% config
ax_cmap={'red','blue'};
str_corr={'$\uparrow\uparrow/\downarrow\downarrow$','$\uparrow\downarrow/\downarrow\uparrow$'};


% vis
for ii=idx_tau_vis
    H=figure('Name','g2_kloc','Units',config_fig.units,'Position',[0 0 7 7],'Renderer',config_fig.rend);
    H.Name=strcat(H.Name,'_',num2str(ii));
    
    ax=tight_subplot(2,1,[0,0],[0.15,0.10],[0.15,0.15]);
    
    axes(ax(1));
    plotFlatMapWrappedRad(gaz_filled,gel_filled,g2_anti{ii},'rect','texturemap');
    
    axes(ax(2));
    plotFlatMapWrappedRad(gaz_filled,gel_filled,g2_corr{ii},'rect','texturemap');
    
    % annotation
    for jj=1:2
        axes(ax(jj));
        ax_pos0 = plotboxpos(ax(jj));
        cbar=colorbar();
        cbar.Position(3)=0.03;
        cbar.Position(1)=sum(ax_pos0([1,3]))+cbar.Position(3);
        cbar.Position([2,4])=ax_pos0([2,4]);
        
        if jj==1
            titlestr = sprintf('%s = %0.3g ms','$\tau$',tau(ii));
            title(titlestr);
            title(cbar,'g^{(2)}');
        end
        
        if jj==1
            colormap(gca,'red');
        elseif jj==2
            colormap(gca,'blue');
        end
        
        set(gca,'Layer','Top');
        axis tight;
        axis equal;
        box on;
        set(gca,'FontSize',config_fig.ax_fontsize);
        set(gca,'LineWidth',config_fig.ax_lwid);
        set(gca,'XTick',-180:90:180);
        set(gca,'YTick',-90:45:90);
        
        if jj==2
            xlabel('$\theta$ (deg)');
        else
            set(gca,'XTickLabel','');
        end
        ylabel('$\phi$ (deg)');
        
        % correlation type
        tt=text(0,67.5,str_corr{jj},'HorizontalAlignment','center','VerticalAlignment','middle','BackgroundColor','w','FontSize',config_fig.ax_fontsize);
        uistack(tt,'bottom');
        
        % hatchfill
        hatch_axis(gca);
    end
end


%% VIS: 1D profiles
for ii=idx_tau_vis
    H=figure('Name','g2_kloc_equator','Units',config_fig.units,'Position',[0 0 7 3],'Renderer',config_fig.rend);
    H.Name=strcat(H.Name,'_',num2str(ii));
    
    hold on;
    
    plot(rad2deg(vaz_filled),g2_anti{ii}(:,idx_el0),'r','DisplayName','$\uparrow\uparrow/\downarrow\downarrow$');
    plot(rad2deg(vaz_filled),g2_corr{ii}(:,idx_el0),'b--','DisplayName','$\uparrow\downarrow$');
    
    % annotate
    box on;
    ax=gca;
    ax.FontSize=config_fig.ax_fontsize;
    ax.LineWidth=config_fig.ax_lwid;
    xlabel('$\theta$ (deg)');
    ylabel('$g^{(2)}$');
    xlim([-180,180]);
    set(gca,'XTick',-180:90:180);
    
%     lgd=legend('Location','best');
end
