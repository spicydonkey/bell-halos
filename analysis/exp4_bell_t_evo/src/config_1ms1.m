% load raw data for efficient data analysis
% correlations under uniform global rotation
%
% 2018.04.28: Time evolution of correlation under pi/2 global rotation
%
% DKS
% 

exp_param=expparams();

%% CONFIGS
% configs.path.base='C:\Users\HE BEC\bell\2018april\exp4_tevo\1.1ms';
% configs.path.data=fullfile(configs.path.base,'data');

configs.path.base='C:\Users\HE BEC\bell_data\exp4_bell_t_evo\1.1ms';
configs.path.data=configs.path.base;

configs.path.out=fullfile(configs.path.base,'out');
configs.path.src=fullfile(configs.path.base,'src');
configs.load.path=fullfile(configs.path.data,'d');

configs.flag.param_scan=0;      % 0 for no scan; 1 for param scan
configs.path.paramlog=fullfile(configs.path.data,'LOG_parameters.txt');
configs.misc.param=1.1;     % DELAY (ms)

configs.load.id=[];
configs.load.mincount=1000;
configs.load.maxcount=Inf;     

configs.load.window{1}=[0.39,0.45];         % T [s]
configs.load.window{2}=[-35e-3,35e-3];      % X [m]
configs.load.window{3}=[-35e-3,35e-3];      % Y [m]


vz=exp_param.vz;

configs.mf(1).mf=1;
configs.mf(1).window={vz*[0.3952,0.4121],35e-3*[-1,1],35e-3*[-1,1]};
configs.mf(1).p_bec0=[1.6229,  -3.1e-3,    -1.5e-3;
                   1.6739,  -3.9e-3,    -1.4e-3];
configs.mf(1).r_bec0=8e-3;

configs.mf(2).mf=0;
configs.mf(2).window={vz*[0.4121,0.4296],35e-3*[-1,1],35e-3*[-1,1]};
configs.mf(2).p_bec0=[1.6877,  0.4e-3,    2.9e-3;
                    1.7420,  0.6e-3,    3.7e-3];
configs.mf(2).r_bec0=8e-3;
                           

% configs.mf{3}.mf=-1;
% configs.mf{3}.bec=[vz*0.4251,  -3.7e-3,    5.1e-3;
%                    vz*0.4381,  -2.0e-3,    5.1e-3];


%%% filter
% 1st stage
configs.filt.r_ball=13e-3;
configs.filt.r_crop=[0.6,1.2];
configs.filt.z_cap=0.8;

% post-distortion
configs.filt2.r_crop=[0.9,1.1];
configs.filt2.z_cap=0.75;

%%% halo centering
configs.post.Dk{1}=[0,0,0];
configs.post.Dk{2}=[0.027,0.007,-0.005];