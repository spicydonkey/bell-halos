function generate_parameter_list()
global pars_fname
%GENERATE_PARAMETER_LIST 
% create PAR_VALUE and write to PARS_FNAME
%creates the list of parameters that will be steped through
%main improvement from v0.2 is the multi parameter feature which stores
%each combination of variables on a new line
%par_value should be fromated
%value1shot1 value2shot1
%value1shot2 value2shot2
%...

%-----------HINT-----------
%Combining variables
%option one 
%product
%a=[1 2]
%b=[5 6]
%transpose(combvec(a,b))
%     1     5
%     2     5
%     1     6
%     2     6
%which maps through all posible combinations

%Option two list
%steps though both in sequenc
%requires both lists to be same length !
%transpose([a; b])
%     1     5
%     2     6
%-----------END OF HINT-----------


% SET FACTORY FUNCTION HERE:
% amp=repelem(linspace(0.05, .35, 50),2);
% duration=linspace(1, 5, 5)*1e-6;
% gmod=[0,1,2,0.5];
% par_value=transpose(combvec(amp,duration,gmod));
%par_value=transpose(amp);
    %par_value=transpose(combvec(amp,duration,gmod));
    %par_value=transpose(combvec(amp,gmod));
    
    %par_value=transpose(linspace(0.738e6-0.2e6,0.738e6+0.2e6,10))
    
%     amp=repmat(linspace(0, 0.25, 400),20);
%    par_value=[[0 1.8e-6 1];[0.1 1.8e-6 1];[0.2 1.8e-6 1]]
%    a=[1 2]*10^6;
%    b=[5 6];
%    par_value=transpose(combvec(a,b,a));

%par_value=transpose(repmat(linspace(0, 0.27, 1600),1,10));

%can then randomize the order of rows, this prevents stripes on plots
%due to long term drift of the exp at the expense of dificulty
%looking at data as it is being collected
%par_value=par_value(randperm(size(par_value,1)),:);

% 2017-10-27: Bell v2, mixing beam characterisation for mf=0
% par_value=linspace(0.1,0.44,10)';

% 2017-10-29: Bell v2
% par_value is a N-by-2 array of [MF, K_R_MIX]
% - spin-pol halos between mf with no rotation pulse to characterise source with noise
% - mf=1 only: more raman amp data
% par_mf_source=[0,0;1,0];
% krmix_list=linspace(0.1,0.44,10);
% par_mf1_rot=combvec(1,krmix_list)';
% 
% par_value=vertcat(par_mf_source,par_mf1_rot);

% 2018-01-22: Bell v2, mixing beam characterisation for Rabi frequency and mF asymmetry
% % par_value is a N-by-2 array of [MF, T_Raman_mix]
% mf=[0,1];
% % traman=1e-6*linspace(0,4,4);
% traman=1e-6*linspace(0,5,6);        % 2018.01.25
% par_value=combvec(mf,traman)';


% 2018-02-19: Bell v2, mixing beam characterisation for Rabi frequency
% % par_value is a N-by-1 array of [T_Raman_mix]
% traman=1e-6*linspace(0,10,11);        % 2018.02.19
% par_value=traman';

% 2018-03.02
% global search
% par_UB_mix=linspace(2.8384e+06,2.9192e+06,50);
% par_value=par_UB_mix';

% global search 4 (~sigma -)
% par_UB_mix=linspace(2.86e+06,2.92e+06,10);        %2.889
% par_T_mix=linspace(10e-6,150e-6,8);
% par_K_mix=linspace(0.2,0.45,8);
% par_Gs_mix=linspace(0,2,3);
% 
% par_value=combvec(par_UB_mix,par_T_mix,par_K_mix,par_Gs_mix)';
% % v4_2: continue from before but "shift" 750 good shots
% n_iter_shift=700;       % we shift by 700 (751 good shots but we can check reproducibility)
% par_value=circshift(par_value,-n_iter_shift);

%%% global - single-pass 
% 2018-03.09
% 1. pi
% par_UB_mix=linspace(1.3e+06,2.1e+06,60);
% par_T_mix=linspace(10e-6,60e-6,5);
% 
% par_value=combvec(par_UB_mix,par_T_mix)';

% % 1. sigma-minus
% par_UB_mix=linspace(1.4e+06,1.5e+06,100);
% par_T_mix=linspace(10e-6,100e-6,10);
% 
% par_value=combvec(par_UB_mix,par_T_mix)';


%%% global - single-pass: XZ-bias B field
% % 2018.03.22
% par_UB_mix=linspace(1.57e+06,1.63e+06,60);
% par_T_mix=linspace(0,20e-6,10);
% par_K_mix=linspace(0.3,0.4,2);
% 
% par_value=combvec(par_UB_mix,par_T_mix,par_K_mix)';
% Irand=randperm(size(par_value,1));      % randomise
% par_value=par_value(Irand,:);


%%% global - prototype Bell
% 2018.03.23
% par_T_mix=1e-6*[1.25,5,10];     % assuming T_pi~10us --> [pi/8,pi/2,pi]
% par_value=par_T_mix';

% % 2018.03.29
% T_pi=10e-6;     % pi-pulse duration
% par_T_mix=T_pi*[1/4, 3/4, 2];     % search phase [pi/4, 3*pi/4, 2*pi]
% par_value=par_T_mix';

%%% IDEAL global rotation
% 2018.04.11 v1
% T_pi=10e-6;     % pi-pulse duration
% par_T_mix=T_pi*[0, 1/2, 1];     % search phase [0, pi/2, pi/4]
% par_value=par_T_mix';

% 2018.04.13 v2
% T_pi=10e-6;     % approx pi-pulse duration
% par_T_mix=T_pi*[1/8,2/8,3/8];     % search phase [pi/8, pi/4, 3*pi/8]
% par_value=par_T_mix';

%------------------------------------------------
% 2018.04.25: triplet global rotation
T_pi=10e-6;     % approx pi-pulse duration
% par_T_mix=T_pi*linspace(0,1,9);     % search phase [0..pi]
% par_value=par_T_mix';

%v3
% short version -- experiment unstable
% par_T_mix=(0.5*T_pi)*[1,2,4,5]*1/6;     % search phase [0..pi/2]
% par_value=par_T_mix';

%v4
% par_T_mix=T_pi*(7:12)/12;     % search phase [pi*7/12..pi]
% par_value=par_T_mix';

% Y-axis 
%   v1
par_T_mix=T_pi*(1:6)/12;     % search phase [0..pi/2]
par_value=par_T_mix';



% %%% char mf=0 rot
% 2018.04.13
% par_T_mix=linspace(0,30e-6,31);
% par_value=par_T_mix';


%Estimate Runtime
cycle_time=25;%in seconds
fprintf(2,'Predicted time to finish %i shots based on cycle time of %2.1f seconds is %s dd HH:MM.\n',size(par_value,1),cycle_time,datestr(datenum(0,0,0,0,0,cycle_time)*size(par_value,1),'dd HH:MM'))

% Save parameter list to file:
%fID = fopen(pars_fname, 'w');
%fprintf(fID,'%d\n', par_value); 
dlmwrite(pars_fname,par_value);
%fclose(fID);
    
end

%% ADD PARAMETER FACTORY FUNCTION HERE:
function K_B=get_K_B()
    N__per_interval=20;
    K_B=[];
    K_B=[K_B, linspace(0, 0.3, N__per_interval)];
    %K_B=[0.1, 0.15, 0.17, 0.19, 0.208];
    
end

function Trigg_Delay=get_Trigg_Delay()
    N__per_val=200;
    Trigg_Delay=[];
    Trigg_Delay=[Trigg_Delay, 1e-3*linspace(0, 20, N__per_val)];    
end

function T_flight=get_T_flight()
    N__per_interval=150;
    T_flight=[];
    T_flight=[T_flight, 1e-6*linspace(0, 2000, N__per_interval)];
    T_flight=[T_flight, 1e-6*linspace(400, 1500, N__per_interval)];
    
%     T_flight=[T_flight, 1e-6*linspace(0, 200, N__per_interval)];
% 
%     T_flight=[T_flight, 1e-6*linspace(200, 500, N__per_interval)];
%     T_flight=[T_flight, 1e-6*linspace(200, 500, N__per_interval)];
% 
%     T_flight=[T_flight, 1e-6*linspace(0, 200, N__per_interval)];
% 
%     T_flight=[T_flight, 1e-6*linspace(600, 2000, N__per_interval)];
%     T_flight=[T_flight, 1e-6*linspace(600, 2000, N__per_interval)];
% 
%     T_flight=[T_flight, 1e-6*linspace(0, 200, N__per_interval)];
% 
%     T_flight=[T_flight, 1e-6*linspace(600, 2000, N__per_interval)];
%     T_flight=[T_flight, 1e-6*linspace(600, 2000, N__per_interval)];
% 
%     T_flight=[T_flight, 1e-6*linspace(0, 200, N__per_interval)];
end

function K_R=get_K_R()
    N__per_interval=40;
    K_R=[];
    
    K_R=[K_R, linspace(0, 1, N__per_interval)];
    K_R=[K_R, linspace(0.5, 0.6, N__per_interval)];
    K_R=[K_R, linspace(0, 1, N__per_interval)];
    K_R=[K_R, linspace(0, 1, N__per_interval)];
end