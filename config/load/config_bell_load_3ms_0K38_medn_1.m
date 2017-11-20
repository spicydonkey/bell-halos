% configures loadHalos

%% Bell experiment configs
configs.exp.name='Bell';
configs.exp.paramNames={'ampRaman','delay','intRaman','comment'};
configs.exp.paramVals={0.38,3,2.5,'medn'}; 
% experimental parameters are: 
%   {Rotation pulse amplitude K_R_mix,...
%   delay between source and rotation pulse (ms),...
%   noise eater intensity setpoint for Raman laser (V),...
%   any other details}
%

%% FLAGS
configs.flags.savedata=1;
configs.flags.savefigs=0;
configs.flags.verbose=1;
configs.flags.graphics=1;

%% MISCELLANEOUS
configs.misc.vel_z=9.8*0.416;    % atom free-fall vert v at detector hit for T-to-Z conversion;
vz=configs.misc.vel_z;

%% FILES
configs.files.path='\\AMPLPC29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\bell_v2\bell\run_3ms_0.38_medn_1';

%% LOAD
configs.load.id=[];               % file id numbers to use for analysis

% file ID and simple pass/fail
configs.load.mincount=800;          % min counts in window - 0 for no min
configs.load.maxcount=Inf;          % max counts in window - Inf for no max

% Detector/trap alignment
configs.load.rot_angle=0.61;

% TXY window - region of interest ( [] --> no crop )
configs.load.window{1}=[0.38,0.44];      % T [s]
configs.load.window{2}=[-35e-3,35e-3];    % X [m]
configs.load.window{3}=[-35e-3,35e-3];    % Y [m]

%% HALO
%%% HALO PARAMS: BEC counts + oscillation removal for broad capture of halos
configs.bec{1}.pos{1}=[1.6524,-3.4e-3,4.6e-3];   % approx condensate locations (z,x,y)
configs.bec{1}.pos{2}=[1.7055,-2.8e-3,4.7e-3];
configs.bec{1}.Rmax=8e-3;      % max condensate sph radius
configs.bec{1}.r_th=8e-3;
configs.halo{1}.dR=[0.2,0.1];   % radial mask fractional width (coarse,fine)
configs.halo{1}.zcap=0.8;       % z-cutoff (fractional wrt radius)
configs.halo{1}.string='$m_F=0$';

configs.bec{2}.pos{1}=[1.5779,-3.3e-3,3.4e-3];   % approx condensate locations (z,x,y)
configs.bec{2}.pos{2}=[1.6283,-2.9e-3,3.8e-3];
configs.bec{2}.Rmax=8e-3;      % max condensate sph radius
configs.bec{2}.r_th=8e-3;
configs.halo{2}.dR=[0.4,0.1];   % radial mask fractional width (coarse,fine)
configs.halo{2}.zcap=0.8;       % z-cutoff (fractional wrt radius)
configs.halo{2}.string='$m_F=1$';

configs.halo{1}.boost=zeros(1,3);
configs.halo{2}.boost=zeros(1,3);
