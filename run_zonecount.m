% Preliminary halo analysis on Bell test experiment
% requires distinguishable-halos repo for halo analysis tools
% 

clear all; clc; close all;

%% set config
% USER SET CONFIG!
config_bell_test;

verbose=configs.flags.verbose;


%% load TXY
txy=load_txy(configs.files.path,configs.load.id,...
            configs.load.window,configs.load.mincount,configs.load.maxcount,...
            configs.load.rot_angle,configs.flags.build_txy,verbose,configs.flags.graphics);
        
zxy=txy2zxy(txy);       % zxy

% TODO - archive the TXY for fast loading

%% capture each halo
% configs
col={'b','r'};

% preallocate
halo_k=cell(2,1);
bec_cent=cell(2,1);
halo_cent=cell(2,1);
ecent=cell(2,1);
erad=cell(2,1);
evecs=cell(2,1);
Nsc=cell(2,1);
dk=cell(2,1);
halo_modeocc=cell(2,1);

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
    this_dz_poles=cellfun(@(c1,c2)abs(c1(1)-c2(1)),this_bec_cent(:,1),this_bec_cent(:,2),'UniformOutput',false);
    this_dz_poles=mean(vertcat(this_dz_poles{:}));    % pole-pole distance (BEC) in Z
    
    this_bool_zcap=cellfun(@(zxy)abs(zxy(:,1))>(0.5*this_dz_poles*configs.halo{ii}.zcap),this_halo_zxy0,'UniformOutput',false);
    
    this_halo_zxy0=cellfun(@(ZXY,BOOL)ZXY(~BOOL,:),this_halo_zxy0,this_bool_zcap,'UniformOutput',false);
    
    % remove zero counts
    % TODO - need to remove these shots from analysis!
    this_bool_empty_halo=(cellfun(@(x)isequal(size(x,1),0),this_halo_zxy0));
    this_halo_zxy0=this_halo_zxy0(~this_bool_empty_halo);
    
    %% Cull aliased hits
    % cull from centred zxy0 halo data
    alias_deadz=100e-9*configs.misc.vel_z;      % dZ in 100 ns
    this_bool_alias=cellfun(@(x) findalias(x,alias_deadz),this_halo_zxy0,'UniformOutput',false);
    this_halo_zxy0_filt=cellfun(@(x,y) x(~y,:),this_halo_zxy0,this_bool_alias,'UniformOutput',false);
    % TODO - problem when halo_zxy0 empty
    
    % TODO - better var management
    this_halo_zxy0_orig=this_halo_zxy0;       % store original
    this_halo_zxy0=this_halo_zxy0_filt;       % alias filtered
    
    % plot
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
    M_rot=this_evecs;   % rotation matrix: X'Y'Z'(principal ellipsoid) --> XYZ
    this_halo_k=cellfun(@(x) (M_rot\x')',this_halo_k,'UniformOutput',false);     % inverse transform
    
    % Stretch to UNIT sphere: unit in collision wave-vector/momenta
    this_halo_k=cellfun(@(x) x./repmat(this_erad',size(x,1),1),this_halo_k,'UniformOutput',false);
    
    % Reverse transform to original/detector axis
    this_halo_k=cellfun(@(x) (M_rot*x')',this_halo_k,'UniformOutput',false);
    
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
    
    %% store this mF results
    halo_k{ii}=this_halo_k;
    
    % misc
    bec_cent{ii}=this_bec_cent;
    halo_cent{ii}=this_halo_cent;
    ecent{ii}=this_ecent;
    erad{ii}=this_erad;
    evecs{ii}=this_evecs;
    
    Nsc{ii}=this_Nsc;
    dk{ii}=this_dk;
    halo_modeocc{ii}=this_halo_modeocc;
end

%% Simple zonal analysis
%%% config zones, binning
nazim=40;
nelev=40;
binwidth=2*sqrt(((2*pi)/nazim)*(pi/nelev));   % TODO - this isn't ideal - may miss counts

%%% build meshgrid of zones
azim=linspace(-pi,pi,nazim);

% scan all elev angle
elev=linspace(-pi/2,pi/2,nelev);

% % scan within z-cap
% elev_lim=asin(configs.halo{1}.zcap);
% elev=linspace(-elev_lim,elev_lim,nelev);

[AZIM,ELEV]=meshgrid(azim,elev);

%%% get counts in zone
% TODO - currently binning with fixed bin width - try Gaussian convolution
nn_halo=cell(2,1);
for ii=1:2
    nn_halo{ii}=halo_zone_density(halo_k{ii},AZIM,ELEV,binwidth);
    nn_halo{ii}=cat(3,nn_halo{ii}{:});  % nazim x nelev x nshot array
end

%%% statistics
% TODO - once the above "bad shot elimination" is fixed, below shouldn't
% need uniformoutput to false
nn_halo_mean=cellfun(@(x)mean(x,3),nn_halo,'UniformOutput',false);
% nn_halo_std=cellfun(@(x)std(x,0,3)/sqrt(size(x,3)),nn_halo);
nn_halo_serr=cellfun(@(x)std(x,0,3)/sqrt(size(x,3)),nn_halo,'UniformOutput',false);

%%% plot
hfig_halo_zone_density=figure(11);
for ii=1:2
    subplot(1,2,ii);
    plot_sph_surf(AZIM,ELEV,nn_halo_mean{ii});
    axis on;
    box on;
    xlabel('$K_X$');
    ylabel('$K_Y$');
    zlabel('$K_Z$');
    title(configs.halo{ii}.string);
    c=colorbar('SouthOutside');
    c.Label.String='Avg. counts';
end