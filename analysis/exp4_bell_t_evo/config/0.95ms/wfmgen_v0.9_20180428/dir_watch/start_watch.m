%Improvements/Bugs
%email notification of eta should inculde days
%ETA calc should be based on the cycle time not the duration since started
%the first two points seem to be the same
%seems that the repeat feature pust the written params out of step with
%what is sent to the waveform gen


%clear all;
global pointer_fname
global DIR_WATCH_ENBL
global watching_dir
global pars_fname
global pars_log_fname
global start_time       %for ETA
global repeats          %number of times the param list has been completed
DIR_WATCH_ENBL=true;

%% Front END:
%START User Variables

%where it the output of the TDC front pannel
watching_dir='Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\';

%set to 1 for Inf no. repeats at end of params
%set to 0 for a single scan
Repeat_Scan=1;

%Low Count Threshold in MiB
%prevents updates if the file is small due to source dropping out
Low_Count_Size=0.05;
%monitoring Email for dropout alert
% Email_add='Bryce.m.Henson@gmail.com';
Email_add='dk.shin1992@gmail.com';
%Email_add='';%use this to disable email alert

%END User Variables
%% Backend

repeats=0;

% create the dir for the pars file and log
%Bryce-i dont think this needs to have a dateappended file for the log
%as all the info should be contained for the log, that way can append
%pars_fname =fullfile(watching_dir,sprintf('%s_parameters.txt', datestr(now,30)));
%pars_log_fname =fullfile(watching_dir,sprintf('%s_LOG_parameters.txt', datestr(now,30)));

% pointer_fname = 'dir_watch/pars/pointer';
pointer_fname =fullfile(watching_dir,sprintf('%s_POINTER.txt', datestr(now,30)));
% pars_fname = 'dir_watch/pars/parameters.txt';
pars_fname =fullfile(watching_dir,sprintf('%s_parameters.txt', datestr(now,30)));
pars_log_fname =fullfile(watching_dir,'LOG_parameters.txt');

%Reset the pointer
fID = fopen(pointer_fname, 'w'); 
fprintf(fID,'%i',1);
fclose(fID);


%create the files
fid3 = fopen(pars_log_fname, 'a');
fclose(fid3);
fid2 = fopen(pars_fname, 'a');
fclose(fid2);

%run
generate_parameter_list();
start_time=clock;
dir_monitor(Low_Count_Size,Email_add,Repeat_Scan)
fclose('all');
