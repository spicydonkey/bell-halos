% Bragg mode occupancy scan
%
%
% 2018.05.19: Bragg pulse scan in trap
%
% DKS
% 

exp_param=expparams();

%% CONFIGS
configs.path.base='C:\Users\David\Documents\bell\bell_data_20180527\exp6_bragg_halo_high_number';

configs.path.data=fullfile(configs.path.base,'data');

configs.path.out=fullfile(configs.path.base,'out');
configs.path.src=fullfile(configs.path.base,'src');

configs.load.path=fullfile(configs.path.data,'d');

configs.flag.param_scan=1;      % 0 for no scan; 1 for param scan
configs.misc.param_const=NaN;      % constant param (no scan)
configs.path.paramlog=fullfile(configs.path.data,'LOG_parameters.txt');


configs.load.id=1:1770;
configs.load.mincount=0;
configs.load.maxcount=Inf;

configs.load.window{1}=[0.37,0.43];         % T [s]
configs.load.window{2}=[-35e-3,35e-3];      % X [m]
configs.load.window{3}=[-35e-3,35e-3];      % Y [m]


vz=exp_param.vz;

configs.mf(1).mf=1;
configs.mf(1).window={vz*[0.3826,0.4001],35e-3*[-1,1],35e-3*[-1,1]};
configs.mf(1).p_bec0=[vz*0.385,   -3.7e-3,    -0.4e-3;
                   vz*0.3969,  -4.4e-3,    -0.3e-3];
configs.mf(1).r_bec0=8e-3;

% configs.mf(2).mf=0;
% configs.mf(2).window={vz*[0.4118,0.429],35e-3*[-1,1],35e-3*[-1,1]};
% configs.mf(2).p_bec0=[1.6861,  0.7e-3,    3.5e-3;
%                     1.7404,  1.0e-3,    4.2e-3];
% configs.mf(2).r_bec0=7e-3;
                           

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


%%% halo centering
configs.post.Dk{1}=[0.009,-0.005,-0.019];

