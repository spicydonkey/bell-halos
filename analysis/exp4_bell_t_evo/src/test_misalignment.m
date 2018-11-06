%% TEST for shape distortion and transform effects on k-resolved correlation
% 
% DKS
% 2018-11-5
%

%% CONFIGS
% data file
%   contains spin,momentum resolved vectors
fdata='C:\Users\HE BEC\Dropbox\phd\data\bell_epr_2018\proc\exp4_tevo\exp4_20181029.mat';

% exp of interest
idx_exp_select=6;          

% vars of interest
vars2load={'k_tau','tau'};

% k-mode of interest
% az_k=pi/2;          
az_k=0;          
el_k=0;
alpha=pi/8;     % cone half-angle

% g2 
% n_dk=7;
% lim_dk=[-0.2,0.2];
n_dk=3;
lim_dk=[-0.1,0.1];

% DISTORTION/TRANSFORM
spin_dK=2;      % SPIN-type to shift (1/2)
dK_dim=3;       % DIM (order: Z,X,Y) for shift
dK_lim=0.08;    
dK_nstep=5;     % ODD to include 0

%% get data
S=load(fdata,vars2load{:});

U=S.k_tau{idx_exp_select};      % vectors (normed k-vecs in zxy coord)
tau=S.tau(idx_exp_select);

clearvars('S');

%% configure shape transorm
% evalute displacements
dK_0=zeros(dK_nstep,3);
dK_0(:,dK_dim)=linspace(-dK_lim,dK_lim,dK_nstep);


%% configure g2
dk_ed_vec=linspace(lim_dk(1),lim_dk(2),n_dk+1);
dk_cent_vec=dk_ed_vec(1:end-1)+0.5*diff(dk_ed_vec);
[~,idx_dk0]=min(abs(dk_cent_vec));    % idx bin nearest zero
dk_ed={dk_ed_vec,dk_ed_vec,dk_ed_vec};
dk_grid_size=cellfun(@(k) length(k)-1,dk_ed);

%% correlations vs. distortion
g2=cell(dK_nstep,1);
g0=cell(dK_nstep,1);
B=NaN(dK_nstep,1);
B0=NaN(dK_nstep,1);

for ii=1:dK_nstep
    % shape transform
    V=U;
    V(:,spin_dK)=boost_zxy(V(:,spin_dK),dK_0(ii,:));
    
    % capture vecs in k-mode
    v_k=cellfun(@(x) inDoubleCone(x,az_k,el_k,alpha),V,'UniformOutput',false);
    
    %% g2
    tg2=summary_disthalo_g2(v_k,dk_ed,0,0,0,0);
    
    tg0=cellfun(@(g) g(idx_dk0,idx_dk0,idx_dk0),tg2);    % BB corr strength
    [tB,tB0]=g2toE(mean(tg0(1:2)),tg0(3));      % correlator
    
    % store
    g2{ii}=tg2;
    g0{ii}=tg0;
    B(ii)=tB;
    B0(ii)=tB0;
end


