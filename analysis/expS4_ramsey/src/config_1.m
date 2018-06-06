% Characterising global rotation around Y-axis
%
%
% DKS
% 2018-05-29


exp_param=expparams();

%% CONFIGS
configs.path.base='C:\Users\HE BEC\bell\2018april\expS4_ramsey';

configs.path.data=fullfile(configs.path.base,'data');

configs.path.out=fullfile(configs.path.base,'out');
configs.path.src=fullfile(configs.path.base,'src');

configs.load.path=fullfile(configs.path.data,'d');

configs.flag.param_scan=1;      % 0 for no scan; 1 for param scan
configs.path.paramlog=fullfile(configs.path.data,'LOG_parameters.txt');

configs.load.id=[];
configs.load.mincount=2500;
configs.load.maxcount=5500;     

configs.load.window{1}=[0.37,0.44];         % T [s]
configs.load.window{2}=[-35e-3,35e-3];      % X [m]
configs.load.window{3}=[-35e-3,35e-3];      % Y [m]


vz=exp_param.vz;

configs.mf(1).mf=1;
configs.mf(1).window={vz*[0.3822,0.3992],35e-3*[-1,1],35e-3*[-1,1]};
configs.mf(1).p_bec0=[vz*0.385,  -2.8e-3,    -1.1e-3;
                   vz*0.3971,  -3.5e-3,    -1.0e-3];
configs.mf(1).r_bec0=6e-3;


configs.mf(2).mf=0;
configs.mf(2).window={vz*[0.3992,0.416],35e-3*[-1,1],35e-3*[-1,1]};
configs.mf(2).p_bec0=[vz*0.4003,  0.9e-3,    3.5e-3;
                    vz*0.4131,  1.1e-3,    4.2e-3];
configs.mf(2).r_bec0=6e-3;
                           

configs.mf(3).mf=-1;
configs.mf(3).window={vz*[0.4166,0.4321],35e-3*[-1,1],35e-3*[-1,1]};
configs.mf(3).p_bec0=[vz*0.4166,  7e-3,    7e-3;
                    vz*0.4302,  7e-3,    7e-3];
configs.mf(3).r_bec0=6e-3;
                           


%%% filter
% 1st stage
configs.filt.r_ball=14e-3;
configs.filt.r_crop=[0.6,1.2];
configs.filt.z_cap=0.8;

% post-distortion
configs.filt2.r_crop=[0.9,1.1];
configs.filt2.z_cap=0.75;


%%% halo centering
configs.post.Dk{1}=[0,0,0];
configs.post.Dk{2}=[0,0,0];
configs.post.Dk{3}=[0,0,0];