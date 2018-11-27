% config for pulse phase delay experiment
%
% 2018.03.26
% 2018.11.27 - revisited
% DKS
% 

exp_param=expparams();
vz=exp_param.vz;


%% CONFIGS
% configs.path.base='\\AMPLPC29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\bell_global\20180326\pulse_phase_delay';
configs.path.base='C:\Users\HE BEC\bell_data\pulse_phase_delay';
configs.path.data=configs.path.base;

configs.path.out=fullfile(configs.path.base,'out');
configs.path.src=fullfile(configs.path.base,'src');
configs.load.path=fullfile(configs.path.data,'d');

configs.flag.param_scan=1;      % 0 for no scan; 1 for param scan
configs.path.paramlog=fullfile(configs.path.base,'LOG_parameters.txt');


configs.load.id=[];
configs.load.mincount=0;
configs.load.maxcount=Inf;


configs.load.window{1}=[0.37,0.44];         % T [s]
configs.load.window{2}=[-35e-3,35e-3];      % X [m]
configs.load.window{3}=[-35e-3,35e-3];      % Y [m]

%%% mf
configs.mf(1).mf=1;
configs.mf(1).window={vz*[0.3807,0.4],35e-3*[-1,1],35e-3*[-1,1]};
configs.mf(1).p_bec0=[1.562,  1.2e-3,    -2.7e-3;
                   1.616,  1.2e-3,    -2.7e-3];
configs.mf(1).r_bec0=13e-3;

configs.mf(2).mf=0;
configs.mf(2).window={vz*[0.4,0.4162],35e-3*[-1,1],35e-3*[-1,1]};
configs.mf(2).p_bec0=[1.638,  -0.3e-3,    2.9e-3;
                    1.69,  -0.3e-3,    2.9e-3];
configs.mf(2).r_bec0=13e-3;


% configs.mf(3).mf=-1;
% configs.mf(3).window={vz*[0.416,0.4349],35e-3*[-1,1],35e-3*[-1,1]};
% configs.mf(3).p_bec0=[1.71,  -0.3e-3,    7.5e-3;
%                     1.763,  0.4e-3,    7.5e-3];
% configs.mf(3).r_bec0=10e-3;


%%% filter
% 1st stage
configs.filt.r_ball=13e-3;
configs.filt.r_crop=[0.6,1.2];
configs.filt.z_cap=0.8;

% post-distortion
configs.filt2.r_crop=[0.9,1.1];
configs.filt2.z_cap=0.8;
