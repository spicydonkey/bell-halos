%% Time evolution of Psi+ halo
% DKS
% 2018-10-25

%% CONFIGS
% data file
fdata='C:\Users\HE BEC\Dropbox\phd\data\bell_epr_2018\proc\exp4_tevo\exp4_20181025.mat';

% capture region
alpha=pi/12;     % cone half-angle
n_az=12;         % azimuthally equispaced between 0 and pi (excl)
n_el=7;         % elev equispaced between -maxphi and maxphi
phi_max=asin(0.75);

% g2
n_dk=15;
lim_dk=[-0.2,0.2];

% bootstrapping
bs_frac=0.2;
bs_nrep=20;

% vis
f_units='normalized';
f_pos=[0.2,0.2,0.2,0.3];
f_ren='painters';

[c,cl,cd]=palette(3);
% c_gray=0.6*ones(1,3);
% line_sty={'-','--',':'};
% mark_typ={'o','s','^'};
% str_ss={'$\vert\!\uparrow\rangle$','$\vert\!\downarrow\rangle$'};
% str_ss={'$m_J = 1$','$m_J = 0$'};
str_ss={'$\uparrow\uparrow$','$\downarrow\downarrow$','$\uparrow\downarrow$'};
mark_siz=7;
line_wid=1.5;
fontsize=12;

%% load data
load(fdata);

%% create halo sections to investigate
Vaz=linspace(0,pi,n_az+1);
Vaz=Vaz(1:end-1);

Vel=linspace(-phi_max,phi_max,n_el);

[vaz,vel]=ndgrid(Vaz,Vel);    % pairs
n_zone=numel(vaz);

%% create g2 params
dk_ed_vec=linspace(lim_dk(1),lim_dk(2),n_dk+1);
dk_cent_vec=dk_ed_vec(1:end-1)+0.5*diff(dk_ed_vec);
[~,idx_dk0]=min(abs(dk_cent_vec));    % idx bin nearest zero
dk_ed={dk_ed_vec,dk_ed_vec,dk_ed_vec};
dk_grid_size=cellfun(@(k) length(k)-1,dk_ed);

%% MAIN
n_tau=length(k_tau);

% preallocate
c4=cell(3,n_tau,n_az,n_el);     % DIM: corr X tau X theta X phi
c3=cell(n_tau,n_az,n_el);       % DIM: tau X theta X phi
m4=NaN(3,n_tau,n_az,n_el);    
m3=NaN(n_tau,n_az,n_el);

g2=c4;
g0=m4;    
B=m3;
B0=m3;

% g2_bs_mu=c4;
% g2_bs_se=c4;
g0_bs_mu=m4;
g0_bs_se=m4;
B_bs_mu=m3;
B0_bs_mu=m3;
B_bs_se=m3;
B0_bs_se=m3;

for ii=1:n_tau
% for ii=2:n_tau      % DEBUG: skip HUGE dataset
    %% 
    k=k_tau{ii};
    n_data_size=size(k,1);
    
    progressbar(0);
    for jj=1:n_zone
        %%
        [iaz,iel]=ind2sub([n_az,n_el],jj);
        taz=vaz(jj);
        tel=vel(jj);
        
        % capture vecs in section
        [~,b_A]=cellfun(@(x) inCone(x,taz,tel,alpha),k,'UniformOutput',false);     % in k
        [~,b_B]=cellfun(@(x) inCone(x,taz+pi,-tel,alpha),k,'UniformOutput',false);    % in -k
        b_AB=cellfun(@(b1,b2) b1|b2,b_A,b_B,'UniformOutput',false);        % in both
        
        k_in=cellfun(@(x,b) x(b,:),k,b_AB,'UniformOutput',false);
        
        %% g2
        tg2=summary_disthalo_g2(k_in,dk_ed,0,0,0,0);
        
        tg0=cellfun(@(g) g(idx_dk0,idx_dk0,idx_dk0),tg2);    % BB corr strength
        [tB,tB0]=g2toE(mean(tg0(1:2)),tg0(3));      % correlator
        
        %%% store
        g2(:,ii,iaz,iel)=tg2(:);
        g0(:,ii,iaz,iel)=tg0;
        B(ii,iaz,iel)=tB;
        B0(ii,iaz,iel)=tB0;
        
        %% BOOTSTRAPPING
        % set up
        bs_nsamp=round(bs_frac*n_data_size);    % # data to sample for bs
        bs_Isamp=cellfun(@(c) randi(n_data_size,[bs_nsamp,1]), cell(bs_nrep,1),...
            'UniformOutput',false);

        %%% run
        % g2
        tg2_bs=cellfun(@(I) summary_disthalo_g2(k_in(I,:),dk_ed,0,0,0,0),bs_Isamp,...
            'UniformOutput',false);
        tg2_bs=cat(1,tg2_bs{:});    % g2 dist from bs
        tg0_bs=cellfun(@(g) g(idx_dk0,idx_dk0,idx_dk0),tg2_bs);    % corr strength
        
        % B
        [tB_bs,tB0_bs]=g2toE(mean(tg0_bs(:,1:2),2),tg0_bs(:,3));
        
        %%% statistics
        % TODO
%         tg2_bs_mu
%         tg2_bs_se

        tg0_bs_mu=mean(tg0_bs,1);
        tg0_bs_se=sqrt(bs_frac)*std(tg0_bs,0,1);
        
        tB_bs_mu=mean(tB_bs,1);
        tB0_bs_mu=mean(tB0_bs,1);
        tB_bs_se=sqrt(bs_frac)*std(tB_bs,0,1);
        tB0_bs_se=sqrt(bs_frac)*std(tB0_bs,0,1);
        
        %%% store
%         g2_bs_mu(:,ii,iaz,iel)=tg2_bs_mu(:);
%         g2_bs_se(:,ii,iaz,iel)=tg2_bs_se(:);
        g0_bs_mu(:,ii,iaz,iel)=tg0_bs_mu;
        g0_bs_se(:,ii,iaz,iel)=tg0_bs_se;
        B_bs_mu(ii,iaz,iel)=tB_bs_mu;
        B0_bs_mu(ii,iaz,iel)=tB0_bs_mu;
        B_bs_se(ii,iaz,iel)=tB_bs_se;
        B0_bs_se(ii,iaz,iel)=tB0_bs_se;
        
        progressbar(jj/n_zone);
    end
end


%% vis: corrected parity - ALL
%%% ALL
[cc,ccl,ccd]=palette(n_zone);

h=figure('Name','triplet_halo_tevo','Units',f_units,'Position',f_pos,'Renderer',f_ren);
hold on;

for ii=1:n_zone
    [iaz,iel]=ind2sub([n_az,n_el],ii);
    
    tp=ploterr(tau,B0(:,iaz,iel),[],B0_bs_se(:,iaz,iel),'-o','hhxy',0);
    set(tp(1),'MarkerSize',mark_siz,'LineWidth',line_wid,...
        'MarkerFaceColor',ccl(ii,:),'Color',cc(ii,:));
    set(tp(2),'LineWidth',line_wid,'Color',cc(ii,:));
    
end

% annotate
box on;
ax=gca;
xlabel('$\tau~[\textrm{ms}]$');
ylabel('Parity $\bar{\mathcal{B}}_{\pi/2}$');
ax.FontSize=fontsize;
ax.LineWidth=1.2;
ylim([-1.2,1.2]);

%% vis: Polar distribution
% config
iaz_disp=[1,4,7];

for ii=1:length(iaz_disp)
    iaz=iaz_disp(ii);       % azim idx to great circle
    figname=sprintf('B0_tevo_polar_%0.0f',rad2deg(Vaz(iaz)));
    
    % figure
    figure('Name',figname,...
        'Units',f_units,'Position',f_pos,'Renderer',f_ren);
    hold on;
    
    pleg=NaN(n_el,1);
    for jj=1:n_el
        tp=ploterr(tau,squeeze(B0(:,iaz,jj)),[],squeeze(B0_bs_se(:,iaz,jj)),'-o','hhxy',0);
        set(tp(1),'MarkerSize',mark_siz,'LineWidth',line_wid,...
            'MarkerFaceColor',ccl(jj,:),'Color',cc(jj,:),'DisplayName',num2str(rad2deg(Vel(jj)),2));
        set(tp(2),'LineWidth',line_wid,'Color',cc(jj,:));
        pleg(jj)=tp(1);
    end
    
    % annotation
    box on;
    ax=gca;
    xlabel('$\tau~[\textrm{ms}]$');
    ylabel('Parity $\bar{\mathcal{B}}_{\pi/2}$');
    ax.FontSize=fontsize;
    ax.LineWidth=1.2;
    ylim([-1.2,1.2]);
    legend(pleg);
    legend(pleg,'Location','EastOutside');
end

%% vis: Equatorial distribution
% config
[~,iel_0]=min(abs(Vel));
el_0=Vel(iel_0);

figname='B0_tevo_eqt';

% figure
figure('Name',figname,...
    'Units',f_units,'Position',f_pos,'Renderer',f_ren);
hold on;

pleg=NaN(n_az,1);
for ii=1:n_az
    tp=ploterr(tau,squeeze(B0(:,ii,iel_0)),[],squeeze(B0_bs_se(:,ii,iel_0)),'-o','hhxy',0);
    set(tp(1),'MarkerSize',mark_siz,'LineWidth',line_wid,...
        'MarkerFaceColor',ccl(ii,:),'Color',cc(ii,:),'DisplayName',num2str(rad2deg(Vaz(ii)),2));
    set(tp(2),'LineWidth',line_wid,'Color',cc(ii,:));
    pleg(ii)=tp(1);
end

% annotation
box on;
ax=gca;
xlabel('$\tau~[\textrm{ms}]$');
ylabel('Parity $\bar{\mathcal{B}}_{\pi/2}$');
ax.FontSize=fontsize;
ax.LineWidth=1.2;
ylim([-1.2,1.2]);
legend(pleg,'Location','EastOutside');