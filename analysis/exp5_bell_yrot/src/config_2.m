% load raw data for efficient data analysis
% correlations under uniform global rotation
%
% 2018.05.04: raw basis-dependent spin correlation - low mode occupancy
%
% DKS
% 

exp_param=expparams();

%% CONFIGS
configs.path.base='C:\Users\HE BEC\bell\2018april\exp5_low_n_theta_v2';

configs.path.data=fullfile(configs.path.base,'data');

configs.path.out=fullfile(configs.path.base,'out');
configs.path.src=fullfile(configs.path.base,'src');

configs.load.path=fullfile(configs.path.data,'d');

configs.flag.param_scan=1;      % 0 for no scan; 1 for param scan
configs.path.paramlog=fullfile(configs.path.data,'LOG_parameters.txt');

configs.load.id=[];
configs.load.mincount=500;
configs.load.maxcount=1500;     

configs.load.window{1}=[0.39,0.45];         % T [s]
configs.load.window{2}=[-35e-3,35e-3];      % X [m]
configs.load.window{3}=[-35e-3,35e-3];      % Y [m]


vz=exp_param.vz;

configs.mf(1).mf=1;
configs.mf(1).window={vz*[0.3952,0.4118],35e-3*[-1,1],35e-3*[-1,1]};
configs.mf(1).p_bec0=[1.6218,  -2.8e-3,    -1.1e-3;
                   1.6730,  -3.5e-3,    -1.0e-3];
configs.mf(1).r_bec0=7e-3;

configs.mf(2).mf=0;
configs.mf(2).window={vz*[0.4118,0.429],35e-3*[-1,1],35e-3*[-1,1]};
configs.mf(2).p_bec0=[1.6867,  0.9e-3,    3.5e-3;
                    1.7407,  1.1e-3,    4.2e-3];
configs.mf(2).r_bec0=7e-3;
                           

% configs.mf{3}.mf=-1;
% configs.mf{3}.bec=[vz*0.4251,  -3.7e-3,    5.1e-3;
%                    vz*0.4381,  -2.0e-3,    5.1e-3];


%%% filter
% 1st stage
configs.filt.r_ball=14e-3;
configs.filt.r_crop=[0.6,1.2];
configs.filt.z_cap=0.8;

% post-distortion
configs.filt2.r_crop=[0.93,1.07];
configs.filt2.z_cap=0.7;

