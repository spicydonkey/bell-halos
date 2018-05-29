% load raw data for efficient data analysis
% correlations under uniform global rotation
%
% 2018.04.13: Experiment with minimal delays between steps
%	* exploring various basis (theta)
%
% DKS
% 

exp_param=expparams();

%% CONFIGS
configs.path.base='C:\Users\HE BEC\bell\2018april\ideal_global_rotation_v3';

configs.path.data=fullfile(configs.path.base,'data');

configs.path.out=fullfile(configs.path.base,'out');
configs.path.src=fullfile(configs.path.base,'src');

configs.load.path=fullfile(configs.path.data,'d');

configs.flag.param_scan=0;      % 0 for no scan; 1 for param scan
configs.path.paramlog=fullfile(configs.path.data,'LOG_parameters.txt');

configs.load.id=280:3118;       % PUSH drifting in first 280 shots
configs.load.mincount=800;
configs.load.maxcount=1600;     

configs.load.window{1}=[0.38,0.46];         % T [s]
configs.load.window{2}=[-35e-3,35e-3];      % X [m]
configs.load.window{3}=[-35e-3,35e-3];      % Y [m]


vz=exp_param.vz;

configs.mf(1).mf=1;
configs.mf(1).window={vz*[0.3922,0.4111],35e-3*[-1,1],35e-3*[-1,1]};
configs.mf(1).p_bec0=[1.6140,  -1.8e-3,    -1.3e-3;
                   1.6649,  -2.6e-3,    -1.1e-3];
configs.mf(1).r_bec0=6e-3;

configs.mf(2).mf=0;
configs.mf(2).window={vz*[0.4111,0.4301],35e-3*[-1,1],35e-3*[-1,1]};
configs.mf(2).p_bec0=[1.6873,  0.8e-3,    2.9e-3;
                    1.7414,  0.8e-3,    3.6e-3];
configs.mf(2).r_bec0=6e-3;
                           

% configs.mf{3}.mf=-1;
% configs.mf{3}.bec=[vz*0.4251,  -3.7e-3,    5.1e-3;
%                    vz*0.4381,  -2.0e-3,    5.1e-3];


%%% filter
% 1st stage
configs.filt.r_ball=10e-3;
configs.filt.r_crop=[0.6,1.2];
configs.filt.z_cap=0.8;

% post-distortion
configs.filt2.r_crop=[0.9,1.1];
configs.filt2.z_cap=0.8;


