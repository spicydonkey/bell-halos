% config for spin-echo experiment in stable B-field
%
% 2018.04.09
% DKS
% 

exp_param=expparams();
vz=exp_param.vz;


%% CONFIGS
% configs.path.base='C:\Users\HE BEC\bell_data\spin_echo_v1';
% configs.path.base='C:\Users\David\Documents\bell\20180329_easter\spin_echo_v1';
configs.path.base='\\AMPLPC29\He BEC Archive\EXPERIMENT-DATA\BELL_MF\20180122_bell-mf-global\bell_global\20180329_easter\spin_echo_v1';
configs.path.paramlog=fullfile(configs.path.base,'LOG_parameters.txt');

configs.load.path=fullfile(configs.path.base,'d');
configs.load.id=[];
configs.load.mincount=0;
configs.load.maxcount=Inf;

configs.load.window={[0.37,0.44],35e-3*[-1,1],35e-3*[-1,1]};

%%% mf
configs.mf(1).mf=1;
configs.mf(1).window={vz*[0.3807,0.4],35e-3*[-1,1],35e-3*[-1,1]};
configs.mf(1).p_bec0=[vz*0.3841,  1.2e-3,    -2.7e-3;
                   vz*0.3971,  1.2e-3,    -2.7e-3];
configs.mf(1).r_bec0=15e-3;

configs.mf(2).mf=0;
configs.mf(2).window={vz*[0.4,0.4165],35e-3*[-1,1],35e-3*[-1,1]};
configs.mf(2).p_bec0=[vz*0.401,  -0.3e-3,    2.9e-3;
                    vz*0.4138,  -0.3e-3,    2.9e-3];
configs.mf(2).r_bec0=6e-3;


% configs.mf(3).mf=-1;
% configs.mf(3).window={vz*[0.416,0.4349],35e-3*[-1,1],35e-3*[-1,1]};
% configs.mf(3).p_bec0=[1.71,  -0.3e-3,    7.5e-3;
%                     1.763,  0.4e-3,    7.5e-3];
% configs.mf(3).r_bec0=10e-3;


%%% filter
% 1st stage
configs.filt.r_ball=10e-3;
configs.filt.r_crop=[0.6,1.2];
configs.filt.z_cap=0.8;

% post-distortion
configs.filt2.r_crop=[0.85,1.1];
configs.filt2.z_cap=0.8;

