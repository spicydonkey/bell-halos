%Categorise mF-species
%

exp_param=expparams();
vz=exp_param.vz;

configs.mf{1}.mf=1;
configs.mf{1}.bec=[vz*0.3868,  -3.3e-3,    3.4e-3;
                   vz*0.3993,  -3.0e-3,    4.3e-3];

configs.mf{2}.mf=0;
configs.mf{2}.bec=[vz*0.4051,  -3.3e-3,    4.6e-3;
                   vz*0.4182,  -3.0e-3,    5.1e-3];
                              
configs.mf{3}.mf=-1;
configs.mf{3}.bec=[vz*0.4251,  -3.7e-3,    5.1e-3;
                   vz*0.4381,  -2.0e-3,    5.1e-3];

configs.r_bec_capt=6e-3;

% cleaning for halo
configs.r_thermal=10e-3;
configs.r_crop=[0.6,1.3];       % min~0.6 ;max~1.3