%% Mode-averaged g2 correlations during time-evolution
% Demo g2 in angular coords
% 
% DKS
% 
% 2018-10-26

%% configs
% FILES
% path_dir='C:\Users\HE BEC\Dropbox\phd\data\bell_epr_2018\proc\exp4_tevo';
path_dir='C:\Users\David\Dropbox\PhD\data\bell_epr_2018\proc\exp4_tevo';
fname_data='exp4_20181025.mat';

% BOOTSTRAPPING
bs_frac=0.2;
bs_nrep=5;

% VIS
f_units='normalized';
% f_pos=[0.2,0.2,0.2,0.3];
f_pos=[0.2,0.2,0.4,0.6];
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
dth=edge2cent(dth_ed);      % angle between k1 and k2 (0 for CL)
thbb=pi-dth;      % angle between -k1 and k2 (0 for BB)
thbb_lim=sort(pi-dth_lim);

%%% g2 analysis
% preallocate
g2=cell(nfiles,3);
G=cell(nfiles,3);
N=cell(nfiles,3);

progressbar(0);
for ii=1:nfiles
%     tk=k_tau{ii};
    
%     tg2=cell(1,3);
%     tG=cell(1,3);
%     tN=cell(1,3);
    
    [g2(ii,:),G(ii,:),N(ii,:)]=g2_ang_disthalo(k_tau{ii},dth_ed,0);
%     [tg2{1},tG{1},tN{1}]=g2_ang(tk(:,1),dth_ed);
%     [tg2{2},tG{2},tN{2}]=g2_ang(tk(:,2),dth_ed);
%     [tg2{3},tG{3},tN{3}]=g2x_ang(tk,dth_ed);
    
%     %store
%     g2(ii,:)=tg2(:);
%     G(ii,:)=tG(:);
%     N(ii,:)=tN(:);
    
    progressbar(ii/nfiles);
end

%% bootstrapping
n_data_size=shotSize(k_tau);

g2_bs=cell(nfiles,1);

progressbar(0);
for ii=1:nfiles
    % set up
    bs_nsamp=round(bs_frac*n_data_size(ii));    % # data to sample for bs
    bs_Isamp=cellfun(@(c) randi(n_data_size(ii),[bs_nsamp,1]), cell(bs_nrep,1),...
        'UniformOutput',false);
    
    %%% run
    % g2
    tk=k_tau{ii};
    
    tg2_bs=cellfun(@(I) g2_ang_disthalo(tk(I,:),dth_ed,0),bs_Isamp,...
        'UniformOutput',false);
    tg2_bs=cellfun(@(g) cat(2,g{:}),tg2_bs,'UniformOutput',false);  % concat corr types along dim2
    tg2_bs=cat(3,tg2_bs{:});        % concat BS along dim3
    
    % store
    g2_bs{ii}=tg2_bs;

    progressbar(ii/nfiles);
end

% statistics
g2_bs_mu=cellfun(@(g) mean(g,3),g2_bs,'UniformOutput',false);
g2_bs_se=cellfun(@(g) sqrt(bs_frac)*std(g,0,3),g2_bs,'UniformOutput',false);

%% gaussian fit to g2
% CONFIGS
% Constraints: MU = 0, OFFSET = 1
fgauss.mdl='y~amp*exp(-1*x^2/(2*sigma^2))+1';
fgauss.cname={'amp','sigma'};
fgauss.fopt=statset('TolFun',10^-10,'TolX',10^-10,'MaxIter',10^6,'UseParallel',0);

sigma0=0.02;

% FIT
% g2_fit=cell(nfiles,3);
g2_fitmdl=cellfun(@(g) fitnlm(thbb,g,fgauss.mdl,[g(1),sigma0],...
    'CoefficientNames',fgauss.cname,'Options',fgauss.fopt),g2,'UniformOutput',false);

% get fit params
fparam=cellfun(@(mdl) mdl.Coefficients.Estimate,g2_fitmdl,'UniformOutput',false);
famp=cellfun(@(x) x(1),fparam);
fsig=cellfun(@(x) x(2),fparam);

% eval fitted model
th_fit=linspace(0,thbb_lim(2),1e3);
g2_fit=cellfun(@(f) feval(f,th_fit),g2_fitmdl,'UniformOutput',false);

%% Correlator - Parity measurement
[B,B0]=g2toE(mean(famp(:,1:2),2),famp(:,3));

%% vis: angular g2
% figures
for ii=1:nfiles
    figname=sprintf('tau_%0.2f',tau(ii));
    
    % draw figure
    figure('Name',figname,'Units',f_units,'Position',f_pos,'Renderer',f_ren);
    hold on;
    p_leg=NaN(3,1);
    for jj=1:3
        % data
%         tp(jj)=plot(th,g2{ii,jj},...
%             'LineStyle','none','LineWidth',line_wid,'Color',c(jj,:),...
%             'MarkerSize',mark_siz,'Marker',mark_typ{jj},'MarkerFaceColor',cl(jj,:),...
%             'DisplayName',str_ss{jj});
        % data w unc
        pp=myploterr(thbb,g2{ii,jj},[],g2_bs_se{ii}(:,jj),mark_typ{jj},c(jj,:));
        set(pp(1),'MarkerSize',mark_siz,'LineWidth',line_wid,'DisplayName',str_ss{jj});
        set(pp(2),'LineWidth',line_wid);
        p_leg(jj)=pp(1);
        
        % fitted
        plot(th_fit,g2_fit{ii,jj},...
            'LineStyle',line_sty{jj},'LineWidth',line_wid,'Color',c(jj,:));
        
        titlestr=sprintf('%s%0.2f ms','$\tau=$',tau(ii));
        title(titlestr);
    end
    
    % annotation
    box on;
    ax=gca;
    ax.FontSize=fontsize;
    ax.LineWidth=ax_lwidth;
    ax.Layer='top';
    xlabel('$\Delta\theta_{\mathrm{BB}}$ (rad)');
    ylabel('$g^{(2)}$');
    lgd=legend(p_leg,'Location','Northeast');
    lgd.FontSize=fontsize-1;
    xlim(thbb_lim);
end