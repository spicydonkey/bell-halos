% Preliminary halo analysis on Bell test experiment
% requires distinguishable-halos repo for halo analysis tools
% 

clear all; clc; close all;

%% CONFIG
override_config=1;
    load_config_default=1;
        config_default_id=1;
    VERBOSE=1;
do_g2_corr=0;


%% load config file
%%% NO MIXING


%%% BELL
% config_bell_run_1;
% config_bell_run_2;
config_bell_run_3;

%%% mF=0
% config_bell_mf_0_mix_0;

%%% mF=1
% config_bell_mf_1_mix_0;
% config_bell_mf_1_mix_7;

verbose=configs.flags.verbose;      % verbose - used in script

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OVERRIDE CONFIG DEFINED FROM FILE
if override_config
    verbose=VERBOSE;
    
    configs.zone.nazim=90;
    configs.zone.nelev=45;
    
    if load_config_default
        run(sprintf('config_default_%d',config_default_id));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% load TXY
[txy,fout]=load_txy(configs.files.path,configs.load.id,...
            configs.load.window,configs.load.mincount,configs.load.maxcount,...
            configs.load.rot_angle,configs.flags.build_txy,verbose,configs.flags.graphics);
       
zxy=txy2zxy(txy);       % zxy

clearvars txy;

% TODO - archive the TXY for fast loading

%% capture each halo
% configs
col={'b','r'};

% preallocate
halo_k=cell(2,1);
halo_zxy0=cell(2,1);
bec_cent=cell(2,1);
halo_cent=cell(2,1);
ecent=cell(2,1);
erad=cell(2,1);
evecs=cell(2,1);
M_rot_ell=cell(2,1);
Nsc=cell(2,1);
dk=cell(2,1);
halo_modeocc=cell(2,1);

% fid_ok=cell(2,1);
bool_empty_halo=cell(2,1);
for ii=1:2
    %% Halo capture
    % capture BECs at poles
    [this_bec_cent,this_bool_bec]=capture_bec(zxy,configs.bec{ii}.pos,configs.bec{ii}.Rmax,0);
    
    % pop bec + thermal counts
    % bool_bec_combined=cellfun(@(x,y)x|y,bool_bec(:,1),bool_bec(:,2),'UniformOutput',false);
    this_bool_bec_combined=cell_horzcat(this_bool_bec);       % horizontal concat the boolean list
    this_bool_bec_combined=cellfun(@(x)sum(x,2)>0,this_bool_bec_combined,'UniformOutput',false);    % row-wise OR
    
    this_halo_zxy=cellfun(@(ZXY,BOOL)ZXY(~BOOL,:),zxy,this_bool_bec_combined,'UniformOutput',false);
    
    % halo centre: midpoint of 'BEC poles'
    this_halo_cent=cellfun(@(c1,c2)(c1+c2)/2,this_bec_cent(:,1),this_bec_cent(:,2),'UniformOutput',false);
    
    % centre to halo
    this_halo_zxy0=cellfun(@(zxy,c0)zxy-repmat(c0,[size(zxy,1),1]),this_halo_zxy,this_halo_cent,'UniformOutput',false);
    
    
    %% radial capture (spherical shell)
    % distance from halo centre
    this_r0=cellfun(@(x)sqrt(sum(x.^2,2)),this_halo_zxy0,'UniformOutput',false);
    this_rlim=configs.halo{ii}.R{1}*(1+configs.halo{ii}.dR{1}*[-1,1]);       % radial limits for sph-shell capture
    
    this_bool_halo=cellfun(@(r)(r<this_rlim(2))&(r>this_rlim(1)),this_r0,'UniformOutput',false);    % determine halo points
        
    % capture halo!
    this_halo_zxy0=cellfun(@(ZXY,BOOL)ZXY(BOOL,:),this_halo_zxy0,this_bool_halo,'UniformOutput',false);
    
    % remove caps
    % TODO - try removing the caps after ellipsoid fit
    this_dz_poles=cellfun(@(c1,c2)abs(c1(1)-c2(1)),this_bec_cent(:,1),this_bec_cent(:,2),'UniformOutput',false);
    this_dz_poles=mean(vertcat(this_dz_poles{:}));    % pole-pole distance (BEC) in Z
    
    this_bool_zcap=cellfun(@(zxy)abs(zxy(:,1))>(0.5*this_dz_poles*configs.halo{ii}.zcap),this_halo_zxy0,'UniformOutput',false);
    
    this_halo_zxy0=cellfun(@(ZXY,BOOL)ZXY(~BOOL,:),this_halo_zxy0,this_bool_zcap,'UniformOutput',false);
    
    %% Filter DLD detector ringing
    dld_deadtime=200e-9;
    
    [this_halo_zxy0,this_bool_ring]=postfilter_dld_ring(this_halo_zxy0,vz*dld_deadtime);
    % ringing summary
    if true     % verbose
        idx_shots_ring=find(cellfun(@(x) sum(x),this_bool_ring));
        warning('Ringing detected in halo: %d/%d shots',numel(idx_shots_ring),numel(this_bool_ring));
    end
    
    %% plot
    figure(9);
    hold on;
    plot_zxy(this_halo_zxy0,1e5,1,col{ii});
    axis equal;
    box on;
    
    %% Elipsoid fit
    efit_flag='';
    [this_ecent,this_erad,this_evecs,this_v,this_echi2]=ellipsoid_fit(circshift(vertcat(this_halo_zxy0{:}),-1,2),efit_flag);
    
    
    %% Tranform to unit sphere (k-space)
    % Initialise k-vector in XYZ coord
    this_halo_k=cellfun(@(x) circshift(x,-1,2),this_halo_zxy0,'UniformOutput',false);
    
    % Centre to fitted ellipsoid
    this_halo_k(:)=boost_zxy(this_halo_k,-this_ecent');
    
    % Transform to ellipsoid principal axis (Rotation)
    this_M_rot=this_evecs;   % rotation matrix: X'Y'Z'(principal ellipsoid) --> XYZ
    this_halo_k=cellfun(@(x) (this_M_rot\x')',this_halo_k,'UniformOutput',false);     % inverse transform
    
    % Stretch to UNIT sphere: unit in collision wave-vector/momenta
    this_halo_k=cellfun(@(x) x./repmat(this_erad',size(x,1),1),this_halo_k,'UniformOutput',false);
    
    % Reverse transform to original/detector axis
    this_halo_k=cellfun(@(x) (this_M_rot*x')',this_halo_k,'UniformOutput',false);
    
    % transform to ZXY system
    this_halo_k=cellfun(@(x) circshift(x,1,2),this_halo_k,'UniformOutput',false);
    
    % plot
    figure(10);
    hold on;
    plot_zxy(this_halo_k,1e5,1,col{ii});
    axis equal;
    box on;
    
    %% Characterise halo
    [this_Nsc,this_dk]=halo_characterise(this_halo_k,configs.halo{ii}.zcap,verbose);
    this_halo_modeocc=halo_mocc(1,0.033,this_Nsc,0.03);
    
    %% Halo centering
    % boost for best g2 BB centering
    this_halo_k=boost_zxy(this_halo_k,configs.halo{ii}.boost);
    
    %% remove shots with zero counts in halo
    % check if any halo is empty
    % NOTE: no more counts should be culled after this stage!
    this_bool_empty_halo=(cellfun(@(x)isequal(size(x,1),0),this_halo_zxy0));
    
    %% store this mF results
    halo_k{ii}=this_halo_k;
    halo_zxy0{ii}=this_halo_zxy0;
    
    % misc
    bool_empty_halo{ii}=this_bool_empty_halo;
    
    bec_cent{ii}=this_bec_cent;
    halo_cent{ii}=this_halo_cent;
    ecent{ii}=this_ecent;
    erad{ii}=this_erad;
    evecs{ii}=this_evecs;
    M_rot_ell{ii}=this_M_rot;    
    
    Nsc{ii}=this_Nsc;
    dk{ii}=this_dk;
    halo_modeocc{ii}=this_halo_modeocc;
end
clearvars zxy;      % clear the huge variable - stores all counts (inc BEC)
clearvars this_*;   % clear all temporary variables

% cull shots with empty halos
bool_empty_halo=or(bool_empty_halo{:});     % combine with OR
fout.id_ok=fout.id_ok(~bool_empty_halo);    % update id_ok
halo_k=cellfun(@(x) x(~bool_empty_halo),halo_k,'UniformOutput',false);

%% Zonal analysis
%%% config zones, binning
nazim=configs.zone.nazim;
nelev=configs.zone.nelev;

binmethod=configs.zone.binmethod;
binwidth=configs.zone.binwidth;

%%% build meshgrid of zones
% TODO - the closed mesh problem
switch binmethod
    case 1
%         azim_vec=linspace(-pi,pi,nazim+1);  % edges
%         
%         % scan over all elev angle
%         elev_vec=linspace(-pi/2,pi/2,nelev+1);
%         % % scan within z-cap
%         % elev_lim=asin(configs.halo{1}.zcap);
%         % elev_vec=linspace(-elev_lim,elev_lim,nelev+1);
%         
%         azim_cent=(azim_vec(1:end-1)+0.5*diff(azim_vec));
%         elev_cent=(elev_vec(1:end-1)+0.5*diff(elev_vec));

        azim_vec=linspace(-pi,pi,nazim);  % edges

        % scan over all elev angle
        elev_vec=linspace(-pi/2,pi/2,nelev);    % easier to set nelev
%         % scan within z-cap
%         elev_lim=asin(configs.halo{1}.zcap);
%         elev_vec=linspace(-elev_lim,elev_lim,nelev);

        azim_cent=azim_vec;
        elev_cent=elev_vec;

    case 2
%         azim_vec=linspace(-pi,pi,nazim+1);
%         azim_vec=azim_vec(2:end);     % center: elminate -pi,pi wrapping
        
        azim_vec=linspace(-pi,pi,nazim);
        
%         % scan over all elev angle
%         elev_cent=linspace(-pi/2,pi/2,nelev);
        % scan within z-cap
        elev_lim=asin(configs.halo{1}.zcap);
        elev_vec=linspace(-elev_lim,elev_lim,nelev);
        
        azim_cent=azim_vec;
        elev_cent=elev_vec;
end
[azim_grid,elev_grid]=meshgrid(azim_cent,elev_cent);

%%% get counts in zone
% TODO - currently binning with fixed bin width - try Gaussian convolution
nn_halo=cell(2,1);
for ii=1:2
    %     nn_halo{ii}=halo_zone_density(halo_k{ii},azim_grid,elev_grid,binwidth,binmethod);
    nn_halo{ii}=halo_zone_density(halo_k{ii},azim_vec,elev_vec,binwidth,binmethod);
    nn_halo{ii}=cat(3,nn_halo{ii}{:});  % nazim x nelev x nshot array
end

%%% statistics
nn_halo_mean=cellfun(@(x)mean(x,3),nn_halo,'UniformOutput',false);
nn_halo_std=cellfun(@(x)std(x,0,3),nn_halo,'UniformOutput',false);
nn_halo_serr=cellfun(@(x)std(x,0,3)/sqrt(size(x,3)),nn_halo,'UniformOutput',false);

%% plot - zone count distribution
hfig_halo_zone_density=figure(11);
for ii=1:2
    subplot(1,2,ii);
    plot_sph_surf(azim_grid,elev_grid,nn_halo_mean{ii});
    axis on;
    box on;
    xlabel('$K_X$');
    ylabel('$K_Y$');
    zlabel('$K_Z$');
    title(configs.halo{ii}.string);
    c=colorbar('SouthOutside');
    c.Label.String='Avg. counts';
    c.TickLabelInterpreter='latex';
    c.Label.Interpreter='latex';
end

%% cull background streak
% mf=1 halo looks fine
max_nn=max(max(nn_halo_mean{2}));   % max counts in zone from mf=1

% background from spontaneous scattering leaves a bright streak in mf=0
nn_thresh=max_nn*2/3;     % TODO - seems to work
bool_thresh=(nn_halo_mean{1}>nn_thresh);    % boolean on meshgrid to treat as noisy

%%% cull mf=0
% set over threshold to NaN
nn_halo{1}(repmat(bool_thresh,[1,1,size(nn_halo{1},3)]))=NaN;

% redo statistics - shouldn't change!
%%% statistics
nn_halo_mean=cellfun(@(x)mean(x,3),nn_halo,'UniformOutput',false);
nn_halo_std=cellfun(@(x)std(x,0,3),nn_halo,'UniformOutput',false);
nn_halo_serr=cellfun(@(x)std(x,0,3)/sqrt(size(x,3)),nn_halo,'UniformOutput',false);

%% plot - background filtered distribution
hfig_halo_zone_density=figure(14);
for ii=1:2
    subplot(1,2,ii);
    plot_sph_surf(azim_grid,elev_grid,nn_halo_mean{ii});
    axis on;
    box on;
    xlabel('$K_X$');
    ylabel('$K_Y$');
    zlabel('$K_Z$');
    title(configs.halo{ii}.string);
    c=colorbar('SouthOutside');
    c.Label.String='Avg. counts';
    c.TickLabelInterpreter='latex';
    c.Label.Interpreter='latex';
end

%% Mixing characterisation
% Rabi cycle coefficient
P_rabi=nn_halo_mean{2}./(nn_halo_mean{1}+nn_halo_mean{2});
% 1 when all state in mf=1; 0 when all states in mf=0

% plot
hfig_rabi_coeff=figure(12);
plot_sph_surf(azim_grid,elev_grid,P_rabi);

axis on;
box on;
xlabel('$K_X$');
ylabel('$K_Y$');
zlabel('$K_Z$');
% title(configs.halo{ii}.string);
c=colorbar('SouthOutside');
c.TickLabelInterpreter='latex';
c.Label.Interpreter='latex';
% c.Label.String='$P(m_F=1)$';
c.Label.String='$P(\uparrow)$';

%% Bell correlations
% TODO - fix all uncertainties in this!

nn_sum=nn_halo{1}+nn_halo{2};   % total number of atoms with same momentum
Jz=nn_halo{2}-nn_halo{1};       % total spin of atoms with same momentum

% simple transform to match back-to-back momentum modes
nn_sum_bb=flip_bb(nn_sum);
Jz_bb=flip_bb(Jz);

% evaluate correlation coefficient
JJ=Jz.*Jz_bb;
NN=nn_sum.*nn_sum_bb;

E=mean(JJ,3)./mean(NN,3);

% TODO
% each grid point has an associated (theta,phi) rotation angle
% Bell inequality in CHSH form is a linear combination of 4 correlation
% coefficients at a precisely related rotation angle configuration.

%% plot correlations
hfig_corr=figure(13);
plot_sph_surf(azim_grid,elev_grid,E);

axis on;
box on;
xlabel('$K_X$');
ylabel('$K_Y$');
zlabel('$K_Z$');
% title(configs.halo{ii}.string);
c=colorbar('SouthOutside');
c.TickLabelInterpreter='latex';
c.Label.Interpreter='latex';
c.Label.String='$E(\theta,\phi)$';

%% calculate g2 correlation
% build halo for g2 corr analysis
nshot=size(halo_k{1},1);
halo_k0=cell(nshot,2);
for ii=1:2
    halo_k0(:,ii)=halo_k{ii};
end
% halo centre correction
for ii=1:2
    halo_k0(:,ii)=boost_zxy(halo_k0(:,ii),configs.halo{ii}.boost);
end

if do_g2_corr
    corr=halo_g2_manager(halo_k0,configs,verbose);
end
