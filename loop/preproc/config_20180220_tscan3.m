%Process raw data
% Summary
%   long pulse - loop scanning mf=1 to get more conclusive data on
%   inefficient Rabi cycle
%
% 2018.02.20
% D K SHIN

configs.load.path='\\AMPLPC29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\d';

configs.load.id=[];             % file id numbers to use for analysis
configs.load.mincount=1500;
configs.load.maxcount=Inf;

configs.load.window{1}=[0.38,0.44];         % T [s]
configs.load.window{2}=[-35e-3,35e-3];      % X [m]
configs.load.window{3}=[-35e-3,35e-3];      % Y [m]