% Configuration file for Bell test - momentum-spin entanglement
%   Characterising single-beam Raman local operation for mF=1 
%
%

%%% FLAGS
configs.flags.do_corr_analysis=0;
configs.flags.force_all_stages=0;    % force all the stages to run (useful for debug)
configs.flags.verbose=2;
configs.flags.savedata=0;       % TODO - req'd currently since each stage passes data by save/load to disk
configs.flags.archive_txy=0;        % archives loaded TXY as .mat file for future reuse
configs.flags.graphics=1;       % toggle to control graphics/plotting options
configs.flags.build_txy=1;

%%% MISCELLANEOUS
configs.misc.vel_z=9.8*0.416;    % atom free-fall vert v at detector hit for T-to-Z conversion;
vz=configs.misc.vel_z;


%% FILES
configs.files.path='\\AMPLPC29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\bell_progress_201709\mix_mf_1\mix_4\d';

% WARNING: MODIFYING BELOW DIR SETTINGS ARE NOT RECOMMENDED
configs.files.dir_data=fileparts(configs.files.path);    % fullpath to data directory
configs.files.archive=fullfile(configs.files.dir_data,'archive');   % dir to archive folder
configs.files.dirout=fullfile(configs.files.dir_data,'output');      % output directory (will be time-stamped)


%% LOAD
configs.load.version=1;         % TXY load stage version number

% file ID and simple pass/fail
configs.load.id=1:436;         % file id numbers to use for analysis
configs.load.mincount=0;         % min counts in window - 0 for no min
configs.load.maxcount=Inf;          % max counts in window - Inf for no max

% Detector/trap alignment
configs.load.rot_angle=0.61;

% TXY window - region of interest ( [] --> no crop )
configs.load.window{1}=[2.442,2.482];      % T [s]
configs.load.window{2}=[-35e-3,35e-3];    % X [m]
configs.load.window{3}=[-35e-3,35e-3];    % Y [m]

%% HALO
%%% HALO PARAMS: BEC counts + oscillation removal for broad capture of halos
configs.bec{1}.pos{1}=[10.0564,-3.4e-3,4.5e-3];   % approx condensate locations (z,x,y)
configs.bec{1}.Rmax{1}=12e-3;      % max condensate sph radius
configs.bec{1}.dR_tail{1}=1;     % BEC tail radial frac diff
configs.bec{1}.pos{2}=[10.1099,-2.9e-3,4.6e-3];
configs.bec{1}.Rmax{2}=12e-3;
configs.bec{1}.dR_tail{2}=1;

configs.halo{1}.R{1}=26e-3;     % estimated radius of halo
configs.halo{1}.dR{1}=0.2;      % broad radial mask fractional width (in/out)
configs.halo{1}.zcap=0.75;   % z-cutoff (fractional wrt radius)
configs.halo{1}.string='$m_F=0$';

configs.bec{2}.pos{1}=[9.9820,-3.7e-3,3.3e-3];   % approx condensate locations (z,x,y)
configs.bec{2}.Rmax{1}=12e-3;      % max condensate sph radius
configs.bec{2}.dR_tail{1}=1;     % BEC tail radial frac diff
configs.bec{2}.pos{2}=[10.0330,-3.0e-3,3.7e-3];
configs.bec{2}.Rmax{2}=12e-3;
configs.bec{2}.dR_tail{2}=1;

configs.halo{2}.R{1}=20e-3;     % estimated radius of halo
configs.halo{2}.dR{1}=0.2;      % broad radial mask fractional width (in/out)
configs.halo{2}.zcap=0.75;   % z-cutoff (fractional wrt radius)
configs.halo{2}.string='$m_F=1$';

configs.halo{1}.boost=zeros(1,3);
configs.halo{2}.boost=zeros(1,3);

%% Spherical zones
configs.zone.nazim=20;
configs.zone.nelev=15;

configs.zone.binmethod=1;
configs.zone.binwidth=2*sqrt(((2*pi)/configs.zone.nazim)*(pi/configs.zone.nelev));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% configs.zone.binmethod
%   1: TODO  simple - divide sphere into simple nonuniform distributed segments
%   (one-to-one). Note: binwidth will not unused.
%   2: angular - requires binwidth (may be one-to-many OR miss counts).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%