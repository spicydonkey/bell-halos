%% Wafeworm generator v0.8
%
%   A simple and easy to use waveform programer for the 33600A 
%       able to produce long constant sections without using up valuable waveform space
%       ability to repeat sub wavefroms
%       built uing sections of code from  Roman Khakimov's 2015 OOP program
%
%   TODO
%       pad wavefrom with a zero
%       determine if it is possible to have sub waveforms with differing sample rates fixed
%
%
%
%--------------------------------------------------------------------------
%                               HINT
%
%%% 1. waveform builder
%   Waveform Format:
%       {'sine',freq(hz),phase(rad),amplitude,Gauss mod,sample rate,duration}
%       {'const',amplitude,sample rate,duration}
%       {'double_sine',freq1(hz),freq2(hz),phase1(rad),phase2(rad),amplitude1,amplitude2,Gauss mod1,Gauss mod2,sample rate,duration}
%
%--------------------------------------------------------------------------


%% CONFIG
doNormalization     =   false;
useVoltageScaling   =   true;
if_plot             =   true;
Trigg_Delay         =   0.0e-3;


% mat_fname='mat_wf/Raman_Bragg.mat';
global max_points srate_max points_min repeats_max
max_points=double(4e6);
srate_max=1e9;              % [S/s] max is 1e9
points_min=double(32);
repeats_max=1e6;
srate_all=srate_max;


f0_AOM=80e6;            % AOM center frequency [Hz]
Ek=84.9e3;              % Bragg diffraction kinetic energy [Hz] 
% NOTE: Ek is completely determined from wavelength (1083.331 nm) and beam geometry (rel angle)
% (plus helium mass and so on)
%   DTheta (deg)    |   Ek
%   90              |   84.9 kHz
%
%


%% Source (first step in Bell)
%%% CONFIG
%%%% Source
% a two-photon Raman pulse with momentum transfer

% XZ-biased trap @0.5G
UB_src=2.549e6;         % Zeeman splitting [Hz] -- in trap - trap bottom
T_src=6.5e-6;             % 7us ~ duration of PI/2 pulse [s]
K_src=0.3;              % max=0.5 (AOM saturation)

Gs_src=3;

phi1=0;
phi2=0;

% evaluate field freqs
dF_Raman_src=-(UB_src-abs(Ek));         % [Hz]    Raman detuning
f1_src=f0_AOM-dF_Raman_src/2;            % [Hz]    Raman L1
f2_src=f0_AOM+dF_Raman_src/2;            % [Hz]    Raman L2


%%% BUILD
%%%% Source (no mixing)
% ch1_params={{'sine', f1_src, phi1, K_src, Gs_src, srate_all, T_src}};
% ch2_params={{'sine', f2_src, phi2, K_src, Gs_src, srate_all, T_src}};


%% Bragg source mF=1 halo
% T_src=16e-6;             % 6us ~ duration of PI/2 pulse [s]
% K_src=0.13;             % max=0.5 (AOM saturation)
% 
% Gs_src=1.8;
% 
% phi1=0;
% phi2=0;
% 
% % evaluate field freqs
% dF_Bragg_src=abs(Ek);         % [Hz]    Bragg
% f1_src=f0_AOM-dF_Bragg_src/2;            % [Hz]    Bragg L1
% f2_src=f0_AOM+dF_Bragg_src/2;            % [Hz]    Bragg L2
% 
% %%% Bragg
% ch1_params={{'sine', f1_src, phi1, K_src, Gs_src, srate_all, T_src}};
% ch2_params={{'sine', f2_src, phi2, K_src, Gs_src, srate_all, T_src}};


%% mF=1 halo + rotation
% %%% Bragg
% T_src=16e-6;             % 6us ~ duration of PI/2 pulse [s]
% K_src=0.13;             % max=0.5 (AOM saturation)
% 
% Gs_src=1.8;
% 
% phi1=0;
% phi2=0;
% 
% % evaluate field freqs
% dF_Bragg_src=abs(Ek);         % [Hz]    Bragg
% f1_src=f0_AOM-dF_Bragg_src/2;            % [Hz]    Bragg L1
% f2_src=f0_AOM+dF_Bragg_src/2;            % [Hz]    Bragg L2
% 
% 

%%% rotation
T_delay_mix=3000e-6;      % Delay between the SRC and MIX pulse
UB_mix=1.599e6;      % Zeeman splitting [Hz] -- at mixing - nuller bias

% T_mix=5e-6;       % SWITCHED
K_mix=0.3;
Gs_mix=1.0;

phi1_mix=pi;
phi2_mix=0;

% evaluate field freqs
dF_Raman_mix=UB_mix;         %[Hz]    Raman detuning
f1_mix=f0_AOM-dF_Raman_mix/2;
f2_mix=f0_AOM+dF_Raman_mix/2;

% exception for 0 duration mix
if T_mix==0
    % cannot create zero-length pulse so we do a 0-amplitude pulse of v-short duration
    K_mix=0;
    T_mix=0.1e-6;
end
% 
% 
% %%% Build
% % Mixing angle characterise (mF=1): Source Bragg + delay + Mixing Raman
% ch1_params={{'sine',f1_src,phi1,K_src,Gs_src,srate_all,T_src},...
%     {'const',0,srate_all,T_delay_mix},...
%     {'double_sine',f1_mix,f2_mix,phi1_mix,phi2_mix,K_mix,K_mix,Gs_mix,Gs_mix,srate_all,T_mix},...
%     };
% 
% ch2_params={{'sine',f2_src,phi2,K_src,Gs_src,srate_all,T_src},...
%     {'const',0,srate_all,T_delay_mix},...
%     {'const',0,srate_all,T_mix},...
%     };
% 
% 
% %%%% REALLY INTERESTING BEHAVIOUR!
% % % TESTING - 2 pulses with phase lag
% % T_bias=1/UB_mix;
% % T_delay_pulse=10*T_bias;
% % phi_delay=pi/4;
% % 
% % ch1_params={{'sine',f1_src,phi1,K_src,Gs_src,srate_all,T_src},...
% %     {'const',0,srate_all,T_delay_mix},...
% %     {'double_sine',f1_mix,f2_mix,phi1_mix,phi2_mix,K_mix,K_mix,Gs_mix,Gs_mix,srate_all,T_mix},...
% %     {'const',0,srate_all,T_delay_pulse},...
% %     {'double_sine',f1_mix,f2_mix,phi1_mix,phi2_mix+phi_delay,K_mix,K_mix,Gs_mix,Gs_mix,srate_all,T_mix},...
% %     };
% % 
% % 
% % ch2_params={{'sine',f2_src,phi2,K_src,Gs_src,srate_all,T_src},...
% %     {'const',0,srate_all,T_delay_mix+2*T_mix+T_delay_pulse},...
% %     };


%% Bell
%%% SOURCE AND ROTATION MUST BE DEFINED

%% BUILD
ch1_params={{'sine', f1_src, phi1, K_src, Gs_src, srate_all, T_src},...
    {'const',0,srate_all,T_delay_mix},...
    {'double_sine',f1_mix,f2_mix,phi1_mix,phi2_mix,K_mix,K_mix,Gs_mix,Gs_mix,srate_all,T_mix},...
    };

ch2_params={{'sine', f2_src, phi2, K_src, Gs_src, srate_all, T_src},...
    {'const',0,srate_all,T_delay_mix},...
    {'const',0,srate_all,T_mix},...
    };



%% Package and send to instrument
%%% Device 1
channels_dev1={ch_to_waveforms(ch1_params),ch_to_waveforms(ch2_params)};
plot_segments(channels_dev1,0);
send_segments(channels_dev1,1);


%%% Device 2 
% TODO 
%   [ ] can't interface with it
%

% channels_dev2={ch_to_waveforms(ch3_raw),ch_to_waveforms(ch4_raw)};
% plot_segments(channels_dev2,0);
% send_segments(channels_dev2,2);    % cant get other function gen running

