% configures run_loop.m

%% FLAGS
configs.flags.savedata=1;
configs.flgas.savefigs=0;
configs.flags.verbose=0;
configs.flags.graphics=1;

%% MISCELLANEOUS
configs.misc.vel_z=9.8*0.416;    % atom free-fall vert v at detector hit for T-to-Z conversion;
vz=configs.misc.vel_z;

%% FILES
configs.files.path_loop='\\amplpc29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\bell_v2\loop';      % local operation path

%% LOAD
% file ID and simple pass/fail
configs.load.mincount=500;         % min counts in window - 0 for no min
configs.load.maxcount=Inf;          % max counts in window - Inf for no max

% Detector/trap alignment
configs.load.rot_angle=0.61;

% TXY window - region of interest ( [] --> no crop )
configs.load.window{1}=[0.38,0.44];      % T [s]
configs.load.window{2}=[-35e-3,35e-3];    % X [m]
configs.load.window{3}=[-35e-3,35e-3];    % Y [m]

%% HALO
% TODO - update with v2 loop info
%%% HALO PARAMS: BEC counts + oscillation removal for broad capture of halos
configs.bec{1}.pos{1}=[vz*0.4049,-3e-3,4.2e-3];   % approx condensate locations (z,x,y)
configs.bec{1}.pos{2}=[vz*0.4178,-2.9e-3,4.6e-3];
configs.bec{1}.Rmax=10e-3;      % max condensate sph radius
configs.bec{1}.dR_tail=1;     % BEC tail radial frac diff
% configs.halo{1}.R{1}=26e-3;     % estimated radius of halo
configs.halo{1}.dR=0.2;      % broad radial mask fractional width (in/out)
configs.halo{1}.zcap=0.8;   % z-cutoff (fractional wrt radius)
configs.halo{1}.string='$m_F=0$';

configs.bec{2}.pos{1}=[vz*0.3861,-3e-3,3.2e-3];   % approx condensate locations (z,x,y)
configs.bec{2}.pos{2}=[vz*0.3988,-3.5e-3,4.2e-3];
configs.bec{2}.Rmax=10e-3;      % max condensate sph radius
configs.bec{2}.dR_tail=1;     % BEC tail radial frac diff
% configs.halo{2}.R{1}=20e-3;     % estimated radius of halo
configs.halo{2}.dR=0.4;      % broad radial mask fractional width (in/out)
configs.halo{2}.zcap=0.8;   % z-cutoff (fractional wrt radius)
configs.halo{2}.string='$m_F=1$';

configs.halo{1}.boost=zeros(1,3);
configs.halo{2}.boost=zeros(1,3);

%% Spherical zones
configs.zone.histtype='gauss';      % {'gauss','latlon'}
configs.zone.nazim=100;
configs.zone.nelev=50;

% configs.zone.histtype='latlon';      % {'gauss','latlon'}
% configs.zone.nazim=100;
% configs.zone.nelev=50;

configs.zone.sig=[0.1,Inf];     % NORM completely integrated
configs.zone.lim=[3,3];

% zone divisions to sample for display
configs.zone.ndiv_az=8;
configs.zone.ndiv_el=8;

% 1D zonal histogram
configs.zone.nbin_dth_1d=300;
configs.zone.g1d_hsize=15;
configs.zone.g1d_sigma=3;