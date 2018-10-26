%% Mode-averaged g2 correlations during time-evolution
% Demo g2 in angular coords
% 
% DKS
% 
% 2018-10-26

%% configs
% FILES
path_dir='C:\Users\HE BEC\Dropbox\phd\data\bell_epr_2018\proc\exp4_tevo';
fname_data='exp4_20181025.mat';

% VIS
f_units='normalized';
f_pos=[0.2,0.2,0.2,0.3];
f_ren='painters';

[c,cl,cd]=palette(3);
% c_gray=0.6*ones(1,3);
line_sty={'-','--',':'};
mark_typ={'o','s','^'};
% str_ss={'$\vert\!\uparrow\rangle$','$\vert\!\downarrow\rangle$'};
% str_ss={'$m_J = 1$','$m_J = 0$'};
str_ss={'$\uparrow\uparrow$','$\downarrow\downarrow$','$\uparrow\downarrow$'};
mark_siz=7;
line_wid=1.5;
fontsize=12;
ax_lwidth=1.2;

%% load data
path_data=fullfile(path_dir,fname_data);
load(path_data);

%% PROTOTYPE: reduce data
warning('Reducing data for speed.');
k_tau{1}=k_tau{1}(1:3000,:);

%% angular g2
% configs
dth_lim=[pi-0.2,pi];
n_dth=30;

% construct diff-angle bins for g2 histogram
dth_ed=linspace(dth_lim(1),dth_lim(2),n_dth);
dth=edge2cent(dth_ed);

%%% g2 analysis
% preallocate
g2=cell(nfiles,3);
G_hist=cell(nfiles,3);
N_hist=cell(nfiles,3);

progressbar(0);
for ii=1:nfiles
    tk=k_tau{ii};
    
    tg2=cell(1,3);
    tG=cell(1,3);
    tN=cell(1,3);
    
    [tg2{1},tG{1},tN{1}]=g2_ang(tk(:,1),dth_ed);
    [tg2{2},tG{2},tN{2}]=g2_ang(tk(:,2),dth_ed);
    [tg2{3},tG{3},tN{3}]=g2x_ang(tk,dth_ed);
    
    %store
    g2(ii,:)=tg2(:);
    G_hist(ii,:)=tG(:);
    N_hist(ii,:)=tN(:);
    
    progressbar(ii/nfiles);
end

%% vis: angular g2
% figures
for ii=1:nfiles
    figname=sprintf('tau_%0.2f',tau(ii));
    
    % draw figure
    figure('Name',figname,'Units',f_units,'Position',f_pos,'Renderer',f_ren);
    hold on;
    tp=NaN(3,1);
    for jj=1:3
        tp(jj)=plot(dth,g2{ii,jj},...
            'LineStyle',line_sty{jj},'LineWidth',line_wid,'Color',c(jj,:),...
            'MarkerSize',mark_siz,'Marker',mark_typ{jj},'MarkerFaceColor',cl(jj,:),...
            'DisplayName',str_ss{jj});
        titlestr=sprintf('%s%0.2f ms','$\tau=$',tau(ii));
        title(titlestr);
    end
    
    % annotation
    box on;
    ax=gca;
    ax.FontSize=fontsize;
    ax.LineWidth=ax_lwidth;
    ax.Layer='top';
    xlabel('$\Delta\theta$ (rad)');
    ylabel('$g^{(2)}$');
    lgd=legend(tp,'Location','Northwest');
    lgd.FontSize=fontsize-1;
    xlim(dth_lim);
end