%% DIAGNOSTIC: pairwise distance histogram
%
%

%% CONFIGS
% data file
%   contains spin,momentum resolved vectors
fdata='C:\Users\HE BEC\Dropbox\phd\data\bell_epr_2018\proc\exp4_tevo\exp4_20181029.mat';

IDX_TAU=4;

% pairwise dist
pair_dist=@(x,y) vnorm(x-y);

dhist.lim=[0,2.3];
dhist.nbin=5e2;
dhist.ed=linspace(dhist.lim(1),dhist.lim(2),dhist.nbin+1);
dhist.cent=edge2cent(dhist.ed);
dhist.binsize=dhist.ed(2)-dhist.ed(1);

% vis
f_units='normalized';
f_pos=[0.2,0.2,0.2,0.3];
f_pos_wide=[0.2,0.2,0.25,0.3];
f_ren='painters';

[c,cl,cd]=palette(3);
% c_gray=0.6*ones(1,3);
line_sty={'-','--',':','-.'};
mark_typ={'o','s','^','d'};
% str_ss={'$\vert\!\uparrow\rangle$','$\vert\!\downarrow\rangle$'};
str_ss={'$m_J = 1$','$m_J = 0$'};
% str_ss={'$\uparrow\uparrow$','$\downarrow\downarrow$','$\uparrow\downarrow$'};
mark_siz=7;
line_wid=1.5;
fontsize=12;
ax_lwidth=1.2;

%% load data
load(fdata);
warning('Reducing data for speed.');
k_tau{1}=k_tau{1}(1:3000,:);


%% main
k=k_tau{IDX_TAU};

% zxy array to vec-sets
V=cell(size(k));
for ii=1:2      % do separately for mJ states
    V(:,ii)=cellfun(@(v) row2set(v), k(:,ii),'UniformOutput',false);
end
% get distances for all pairwise
D=cellfun(@(vc) cell2mat(pairfun(pair_dist, vc)),V,'UniformOutput',false);

% histogram - unique (unordered) and proper (not-self) pairs
ND=cellfun(@(d) nhist( d( triu( true(size(d)) ,1)), {dhist.ed}),...
    D,'UniformOutput',false);

% tidy
n_d=cell(1,2);
for ii=1:2
    n_d{ii}=cat(2,ND{:,ii});
end

%% statistics
n_samp=size(k,1);
n_d_avg=cellfun(@(x) mean(x,2),n_d,'UniformOutput',false);
n_d_se=cellfun(@(x) (1/sqrt(n_samp))*std(x,[],2),n_d,'UniformOutput',false);

n_det_avg=mean(shotSize(k),1);

%% vis
h=figure('Name','pair_dists','Units',f_units,'Position',f_pos,'Renderer',f_ren);
hold on;

H={};
for ii=1:2 
    H{ii}=shadedErrorBar(dhist.cent,n_d_avg{ii},n_d_se{ii},...
        {'Color',c(ii,:),'LineWidth',line_wid});
    H{ii}.mainLine.DisplayName=str_ss{ii};
end
H=[H{:}];

% summary of params
titlestr=sprintf('%s%0.2f ms, %s(%0.2g, %0.2g), %s%d',...
    '$\tau=$',tau(IDX_TAU),...
    '$\bar{N}=$',n_det_avg,...
    '$n_{\mathrm{bin}}=$',dhist.nbin);
title(titlestr);

% annotate
ax=gca;
set(ax,'Layer','Top');
ax.FontSize=fontsize;
box on;
xlabel('$\Delta k$');
ylabel('Avg \# pairs in bin');
lgd=legend([H.mainLine],'Location','NorthWest');
lgd.FontSize=fontsize;