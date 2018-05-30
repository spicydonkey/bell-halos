% (c) Roman Khakimov, 2015
function instr_update(param)
% INSTR_UPDATE(PARAM)
%
% assign WaveformGenMain waveform parameters from input PARAM
% NOTE: this function must be configured by user in each scenario
% Remember to remove variable declaration in WaveformGenMain.m
%

global DIR_WATCH_ENBL

%     T_flight=param
%     K_R=param
%     Trigg_Delay=param
    %B_trap_bottom=param(1)
%     K_R=param(1)
%     %K_R=param(1)
%     T_Raman=param(2)
%     Gs_mod_R=param(3)

%% 2017-10-27: Bell v2, mixing beam characterisation for mf=0
% K_R_mix=param;

%% 2017-10-29: Bell v2
% % par_value is a N-by-2 array of [MF, K_R_MIX]
% % - spin-pol halos between mf with no rotation pulse to characterise source with noise
% % - mf=1 only: more raman amp data
% source_mf=param(1);
% K_R_mix=param(2);

%% 2018-01-22: Bell v2, mixing beam characterisation for Rabi frequency and mF asymmetry
% % par_value is a N-by-2 array of [MF, T_Raman_mix]
% source_mf=param(1);
% T_Raman_mix=param(2);


%% 2018-02-19: Bell v2, mixing beam characterisation for Rabi frequency
% par_value is a N-by-2 array of [MF, T_Raman_mix]
% T_Raman_mix=param(1);


%% 2018-03-02
% UB_mix=param(1);

% 4
% par_value=combvec(par_UB_mix,par_T_mix,par_K_mix,par_Gs_mix)';
% UB_mix=param(1);
% T_mix=param(2);
% K_mix=param(3);
% Gs_mix=param(4);


%% 2018-03-09
% par_value=combvec(par_UB_mix,par_T_mix,par_K_mix,par_Gs_mix)';
% UB_mix=param(1);
% T_mix=param(2);


%% global - single-pass: XZ-bias B field
% 20180323-201804...
T_mix=param(1);



%% generate waveform
%     WaveformGenMainv4;
WaveformGenMain;
end