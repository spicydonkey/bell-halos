% Configuration file for Bell test - momentum-spin entanglement
%
%
clearvars configs;

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
configs.misc.deadtime=300e-9;

%% FILES
configs.files.path='\\AMPLPC29\He BEC Archive\EXPERIMENT-DATA\xstate_mom_corr\90deg_raman_beams\6_mednum\d';

% WARNING: MODIFYING BELOW DIR SETTINGS ARE NOT RECOMMENDED
configs.files.dir_data=fileparts(configs.files.path);    % fullpath to data directory
configs.files.archive=fullfile(configs.files.dir_data,'archive');   % dir to archive folder
configs.files.dirout=fullfile(configs.files.dir_data,'output');      % output directory (will be time-stamped)


%% LOAD
configs.load.version=3;         % TXY load stage version number

% file ID and simple pass/fail
configs.load.id=1:2848;         % file id numbers to use for analysis
configs.load.mincount=1400;         % min counts in window - 0 for no min
configs.load.maxcount=1800;          % max counts in window - Inf for no max

% Detector/trap alignment
configs.load.rot_angle=0.61;

% TXY window - region of interest ( [] --> no crop )
configs.load.window{1}=[2.44,2.48];      % T [s]
configs.load.window{2}=[-35e-3,35e-3];    % X [m]
configs.load.window{3}=[-35e-3,35e-3];    % Y [m]

%% HALO
% configs.bec.pos{1}=[10.0535,-3e-3,4.2e-3];   % approx condensate locations (z,x,y)
% configs.bec.Rmax{1}=6e-3;      % max condensate sph radius
% configs.bec.dR_tail{1}=1.75;     % BEC tail radial frac diff
% configs.bec.pos{2}=[10.0227,-2.1e-3,0.6e-3];
% configs.bec.Rmax{2}=6e-3;
% configs.bec.dR_tail{2}=1.4;
% 
% configs.halo.R{1}=26e-3;     % estimated radius of halo
% configs.halo.dR{1}=0.2;      % broad radial mask fractional width (in/out)
% configs.halo.R{2}=24e-3;
% configs.halo.dR{2}=0.2;
% 
% configs.halo.zcap=0.8;   % z-cutoff (fractional wrt radius)
% 
% % TODO - does boost need to be optimised for different g2 analysis?
% %   currently SINGLE boost applied to halo2 to obtain best g2_01_BB
% configs.halo.boost{1}=zeros(1,3);
% configs.halo.boost{2}=[0.05,0.0,0.00];

configs.bec{1}.pos{1}=[1.00560e+01,-3.3e-3,4.2e-3];   % approx condensate locations (z,x,y)
configs.bec{1}.Rmax{1}=6e-3;      % max condensate sph radius
configs.bec{1}.dR_tail{1}=1.5;     % BEC tail radial frac diff

configs.halo{1}.R{1}=26e-3;     % estimated radius of halo
configs.halo{1}.dR{1}=0.25;      % broad radial mask fractional width (in/out)
configs.halo{1}.zcap=0.8;   % z-cutoff (fractional wrt radius)
configs.halo{1}.string='$m_F=0$';

configs.bec{2}.pos{1}=[1.00250e+01,-2.3e-3,0.6e-3];   % approx condensate locations (z,x,y)
configs.bec{2}.Rmax{1}=6e-3;      % max condensate sph radius
configs.bec{2}.dR_tail{1}=1;     % BEC tail radial frac diff

configs.halo{2}.R{1}=24e-3;     % estimated radius of halo
configs.halo{2}.dR{1}=0.3;      % broad radial mask fractional width (in/out)
configs.halo{2}.zcap=0.8;   % z-cutoff (fractional wrt radius)
configs.halo{2}.string='$m_F=1$';

configs.halo{1}.boost=zeros(1,3);
configs.halo{2}.boost=[0.05,0.0,0.00];


%% Spherical zones
configs.zone.nazim=90;
% zcap=configs.halo.zcap;
% configs.zone.nelev=round((asin(zcap)/(pi/2))*configs.zone.nazim/2);
configs.zone.nelev=45;

configs.zone.binmethod=1;
configs.zone.binwidth=2*sqrt(((2*pi)/configs.zone.nazim)*(pi/configs.zone.nelev));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% configs.zone.binmethod
%   1: TODO  simple - divide sphere into simple nonuniform distributed segments
%   (one-to-one). Note: binwidth will not unused.
%   2: angular - requires binwidth (may be one-to-many OR miss counts).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Bell 


%% g2 corr


% 1) X-halo Cart BB
configs.corr{1}.type.comp=[1,2];           % components to analysis: cross halo 1,2
configs.corr{1}.type.coord='cart';         % Cartesian (ZXY)
configs.corr{1}.type.opt='BB';             % BB / CL
configs.corr{1}.lim=0.2*repmat([-1,1],[3,1]);
configs.corr{1}.nBin=15*[1,1,1];   % number of bins

% % 2) Single-halo cart BB - m_J=0
% configs.corr{2}.type.comp=1;
% configs.corr{2}.type.coord='cart';
% configs.corr{2}.type.opt='BB';
% configs.corr{2}.lim=0.2*repmat([-1,1],[3,1]);
% configs.corr{2}.nBin=11*[1,1,1];   % number of bins
% 
% % 3) Single-halo cart BB - m_J=1
% configs.corr{3}.type.comp=2;
% configs.corr{3}.type.coord='cart';
% configs.corr{3}.type.opt='BB';
% configs.corr{3}.lim=0.2*repmat([-1,1],[3,1]);
% configs.corr{3}.nBin=11*[1,1,1];   % number of bins