%% TESTING: automated centering of BB regions
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
n_dk=15;
lim_dk=[-0.2,0.2];


%% get data
S=load(fdata,vars2load{:});

U=S.k_tau{idx_exp_select};      % vectors (normed k-vecs in zxy coord)
tau=S.tau(idx_exp_select);

clearvars('S');

% get k-mode
V_k=cellfun(@(x) inDoubleCone(x,az_k,el_k,alpha),U,'UniformOutput',false);

%% configure g2
dk_ed_vec=linspace(lim_dk(1),lim_dk(2),n_dk+1);
dk_cent_vec=dk_ed_vec(1:end-1)+0.5*diff(dk_ed_vec);
[~,idx_dk0]=min(abs(dk_cent_vec));    % idx bin nearest zero
dk_ed={dk_ed_vec,dk_ed_vec,dk_ed_vec};
dk_grid_size=cellfun(@(k) length(k)-1,dk_ed);

%% main
% 
% for ++/--
%   eval g2BB(DK)
%   find DK : max(g2BB(DK))
%   centralise by shifting k-vecs by -DK/2
%
% for +-
%   TODO
%

V_k0=V_k;
dk_max={[0,0,0],[0,0,0]};   % initialise
dK=cellfun(@(x) -x/2,dk_max,'UniformOutput',false);

while true
    V_k0=k_recenter(V_k0,dK);
    g2=summary_disthalo_g2(V_k0,dk_ed,0,1,0,0);
    
    [g2max,Imax]=cellfun(@(g) (max(g(:))),g2,'UniformOutput',false);
    
    g2max{:}
    
    sub=cell(1,2);
    for ii=1:2
        sub{ii}=cell(1,ndims(g2{ii}));
        [sub{ii}{:}]=ind2sub(size(g2{ii}),Imax{ii});
    end
    
    dk_max=cellfun(@(s) cellfun(@(I) dk_cent_vec(I),s),sub,'UniformOutput',false);
    dK=cellfun(@(x) -x/2,dk_max,'UniformOutput',false);
    
    dK{:}
    
    drawnow;
end
