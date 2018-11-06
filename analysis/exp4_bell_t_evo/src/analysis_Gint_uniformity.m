%% ANALYSIS: integrated g2(dk) (cartesian) and "uniform" alignment
% DKS
% 2018-11-06

%% CONFIGS
% data file
%   contains spin,momentum resolved vectors
fdata='C:\Users\HE BEC\Dropbox\phd\data\bell_epr_2018\proc\exp4_tevo\exp4_20181029.mat';

% k-modes
alpha=pi/8;                 % cone half-angle
lim_az=[0,pi];              % limit of azim angle (exc +pi)
phi_max=pi/4;
lim_el=[-phi_max,phi_max];
n_az=24;                	% equispaced bins
n_el=15;
az_disp=deg2rad([0,45,90]);     % azim sections (great circles) to display

% g2 
n_dk=15;
lim_dk=[-0.2,0.2];

% integration
dk_int_lim=0.1;

% bootstrapping
bs_frac=0.2;
bs_nrep=20;

%% configure g2
dk_ed_vec=linspace(lim_dk(1),lim_dk(2),n_dk+1);
dk_cent_vec=dk_ed_vec(1:end-1)+0.5*diff(dk_ed_vec);
[~,idx_dk0]=min(abs(dk_cent_vec));    % idx bin nearest zero
dk_ed={dk_ed_vec,dk_ed_vec,dk_ed_vec};
dk_grid_size=cellfun(@(k) length(k)-1,dk_ed);

% integrable region config
dIk_int=round(dk_int_lim/(range(lim_dk)/n_dk));
intFilter=@(x) x(idx_dk0-dIk_int:idx_dk0+dIk_int,...
    idx_dk0-dIk_int:idx_dk0+dIk_int,idx_dk0-dIk_int:idx_dk0+dIk_int);


%% load data
load(fdata);

%% PROTOTYPE: reduce data
warning('Reducing data for speed.');
k_tau{1}=k_tau{1}(1:3000,:);

% warning('DEBUG ONLY!!!!');
% k_tau=cellfun(@(k) k(1:100,:),k_tau,'UniformOutput',false);

%% k-modes
Vaz=linspace(lim_az(1),lim_az(2),n_az+1);
Vaz=Vaz(1:end-1);       % exclude last
Vel=linspace(lim_el(1),lim_el(2),n_el);

[vaz,vel]=ndgrid(Vaz,Vel);    % AZ-EL grid
n_zone=numel(vaz);

% special angles
[~,iaz_disp]=arrayfun(@(x) min(abs(Vaz-x)),az_disp);     % idx to ~displayable azim angles
[~,iel_0]=min(abs(Vel));        % idx to ~zero elev angle (equator)


%% main
n_tau=length(k_tau);

G_int=NaN(n_tau,n_az,n_el,3);
G_sum=NaN(n_tau,n_az,n_el);
B=NaN(n_tau,n_az,n_el);
B0=NaN(n_tau,n_az,n_el);

G_int_se=NaN(n_tau,n_az,n_el,3);
G_sum_se=NaN(n_tau,n_az,n_el);
B_se=NaN(n_tau,n_az,n_el);
B0_se=NaN(n_tau,n_az,n_el);

for ii=1:n_tau
    k=k_tau{ii};
    n_data_size=size(k,1);
    
    progressbar(0);
    for jj=1:n_zone
        %%
        [iaz,iel]=ind2sub([n_az,n_el],jj);
        taz=vaz(jj);
        tel=vel(jj);
        
        k_cone=cellfun(@(x) inDoubleCone(x,taz,tel,alpha),k,'UniformOutput',false);
        
        %% g2
        tg2=summary_disthalo_g2(k_cone,dk_ed,0,0,0,0);
        
        % volume integrated G2
        tG_int=cellfun(@(g) meanall(intFilter(g),'omitnan'),tg2);
        tG_sum=0.5*sum(tG_int([1 2 3 3])-1);
        [tB,tB0]=g_to_E(tG_int);
    
        %% BOOTSTRAPPING
        % set up
        bs_nsamp=round(bs_frac*n_data_size);    % # data to sample for bs
        bs_Isamp=cellfun(@(c) randi(n_data_size,[bs_nsamp,1]), cell(bs_nrep,1),...
            'UniformOutput',false);

        %%% run
        % g2
        tg2_bs=cellfun(@(I) summary_disthalo_g2(k_cone(I,:),dk_ed,0,0,0,0),bs_Isamp,...
            'UniformOutput',false);
        tg2_bs=cat(1,tg2_bs{:});    % g2 dist from bs
        
        % volume integrated G2
        tG_int_bs=cellfun(@(g) meanall(intFilter(g),'omitnan'),tg2_bs);
        tG_sum_bs=0.5*sum([tG_int_bs,tG_int_bs(:,3)],2)-2;
        [tB_bs,tB0_bs]=g_to_E(tG_int_bs);
        
        % stat
        tG_int_mu=mean(tG_int_bs,1);
        tG_int_se=sqrt(bs_frac)*std(tG_int_bs,[],1);
        tG_sum_mu=mean(tG_sum_bs,1);
        tG_sum_se=sqrt(bs_frac)*std(tG_sum_bs,[],1);
        tB_mu=mean(tB_bs,1);
        tB_se=sqrt(bs_frac)*std(tB_bs,[],1);
        tB0_mu=mean(tB0_bs,1);
        tB0_se=sqrt(bs_frac)*std(tB0_bs,[],1);
        
        
        %% store
        G_int(ii,iaz,iel,:)=tG_int;
        G_sum(ii,iaz,iel)=tG_sum;
        B(ii,iaz,iel)=tB;
        B0(ii,iaz,iel)=tB0;
        
        G_int_se(ii,iaz,iel,:)=tG_int_se;
        G_sum_se(ii,iaz,iel)=tG_sum_se;
        B_se(ii,iaz,iel)=tB_se;
        B0_se(ii,iaz,iel)=tB0_se;
        
        progressbar(jj/n_zone);
    end
end