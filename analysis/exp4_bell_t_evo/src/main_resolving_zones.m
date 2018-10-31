%% Halo averaged correlator vs. Resolving momentum modes
% DKS
% 20181031

%% configs
% FILES
path_dir='C:\Users\HE BEC\Dropbox\phd\data\bell_epr_2018\proc\exp4_tevo\';
fname_data='out_20181031.mat';

%%% analysis
% to be loaded from data --> identical to k-resolved analysis
% g2
n_dk=15;         % 15; 29;
lim_dk=[-0.2,0.2];

% bootstrapping
bs_frac=0.2;
bs_nrep=5;      % 20

% vis
f_units='normalized';
f_pos=[0.2,0.2,0.2,0.3];
% f_pos=[0.2,0.2,0.4,0.6];
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

%% load proc data
path_data=fullfile(path_dir,fname_data);
S=load(path_data);

% get data
tau=S.tau;
k_tau=S.k_tau;

%% main
%%% create g2 params
dk_ed_vec=linspace(lim_dk(1),lim_dk(2),n_dk+1);
dk_cent_vec=dk_ed_vec(1:end-1)+0.5*diff(dk_ed_vec);
[~,idx_dk0]=min(abs(dk_cent_vec));    % idx bin nearest zero
dk_ed={dk_ed_vec,dk_ed_vec,dk_ed_vec};
dk_grid_size=cellfun(@(k) length(k)-1,dk_ed);

%%% corr analysis
n_tau=length(k_tau);

C=cell(n_tau,3);
M3=NaN(n_tau,3);
M1=NaN(n_tau,1);

g2=C;
g0=M3;
B=M1;
B0=M1;

g2_bs_mu=C;
g2_bs_se=C;
g0_bs_mu=M3;
g0_bs_se=M3;
B_bs_mu=M1;
B0_bs_mu=M1;
B_bs_se=M1;
B0_bs_se=M1;

progressbar(0);
for ii=1:n_tau
    k=k_tau{ii};
    n_data_size=size(k,1);
    
    %% g2 - halo averaged
    tg2=summary_disthalo_g2(k,dk_ed,0,0,0,0);
    
    tg0=cellfun(@(g) g(idx_dk0,idx_dk0,idx_dk0),tg2);    % BB corr strength
    [tB,tB0]=g2toE(mean(tg0(1:2)),tg0(3));              % correlator
    
    %% bootstrapping
    % set up
    bs_nsamp=round(bs_frac*n_data_size);    % # data to sample for bs
    bs_Isamp=cellfun(@(c) randi(n_data_size,[bs_nsamp,1]), cell(bs_nrep,1),...
        'UniformOutput',false);
    
    % run
    tg2_bs=cellfun(@(I) summary_disthalo_g2(k(I,:),dk_ed,0,0,0,0),bs_Isamp,...
        'UniformOutput',false);
    tg2_bs=cat(1,tg2_bs{:});    % g2 dist from bs
    tg0_bs=cellfun(@(g) g(idx_dk0,idx_dk0,idx_dk0),tg2_bs);    % corr strength
    [tB_bs,tB0_bs]=g2toE(mean(tg0_bs(:,1:2),2),tg0_bs(:,3));
    
    % statistics
    tg2_bs_cat=cell(1,3);       % tidy bs sampled g2
    for kk=1:3
        tg2_bs_cat{kk}=cat(4,tg2_bs{:,kk});
    end
    tg2_bs_mu=cellfun(@(x) mean(x,4),tg2_bs_cat,'UniformOutput',false);
    tg2_bs_se=cellfun(@(x) sqrt(bs_frac)*std(x,0,4),tg2_bs_cat,'UniformOutput',false);
    
    tg0_bs_mu=mean(tg0_bs,1);
    tg0_bs_se=sqrt(bs_frac)*std(tg0_bs,0,1);
    
    tB_bs_mu=mean(tB_bs,1);
    tB0_bs_mu=mean(tB0_bs,1);
    tB_bs_se=sqrt(bs_frac)*std(tB_bs,0,1);
    tB0_bs_se=sqrt(bs_frac)*std(tB0_bs,0,1);
    
    %% store
    g2(ii,:)=tg2;
    g0(ii,:)=tg0;
    B(ii)=tB;
    B0(ii)=tB0;
    
    g2_bs_mu(ii,:)=tg2_bs_mu;
    g2_bs_se(ii,:)=tg2_bs_se;
    g0_bs_mu(ii,:)=tg0_bs_mu;
    g0_bs_se(ii,:)=tg0_bs_se;
    
    B_bs_mu(ii)=tB_bs_mu;
    B0_bs_mu(ii)=tB0_bs_mu;
    B_bs_se(ii)=tB_bs_se;
    B0_bs_se(ii)=tB0_bs_se;
    
    progressbar(ii/n_tau);
end

%% vis: Parity
% config
az_disp=deg2rad([0,45,90]);

% FIGURE
h=figure('Name','parity_tevo','Units',f_units,'Position',f_pos,'Renderer',f_ren);
hold on;

% VOLUME AVERAGED
tp=myploterr(tau,B0,[],B0_bs_se,'d',[0,0,0]);
set(tp(1),'MarkerSize',mark_siz,'LineWidth',line_wid,'DisplayName','$V$ averaged');
set(tp(2),'LineWidth',line_wid);
pleg_avg=tp(1);

% MODE RESOLVED (equator only)
pleg_res=NaN(numel(az_disp),1);
[~,iaz_disp]=arrayfun(@(x) min(abs(S.Vaz-x)),az_disp);     % idx to ~displayable azim angles
[~,iel_0]=min(abs(S.Vel));        % idx to ~zero elev angle (equator)

for ii=1:numel(az_disp)
    iaz=iaz_disp(ii);
    taz=S.Vaz(iaz);
    pname=sprintf('%s %0.2g','$\theta=$',rad2deg(S.Vaz(iaz_disp(ii))));
    
    tp=myploterr(tau,S.B0(:,iaz,iel_0),[],S.B0_bs_se(:,iaz,iel_0),'o',c(ii,:));
    set(tp(1),'MarkerSize',mark_siz,'LineWidth',line_wid,...
        'DisplayName',pname);
    set(tp(2),'LineWidth',line_wid);
    pleg_res(ii)=tp(1);
end

% annotation
box on;
ax=gca;
set(ax,'Layer','Top');
xlabel('$\tau~[\textrm{ms}]$');
ylabel('Parity $\bar{\mathcal{B}}_{\pi/2}$');
ax.FontSize=fontsize;
ax.LineWidth=1.2;
ylim([-1.2,1.2]);
lgd=legend([pleg_avg;pleg_res],'Location','SouthWest');
lgd.FontSize=fontsize-1;
