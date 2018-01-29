%Process raw data
% Summary
%   First attempt at investigating mF asymmetry in Raman rotation 
%   loop data scanning 
%
% 2018.01.29
% D K SHIN

configs.load.path='C:\Users\HE BEC\exp-data\bell\loop_tscan\d';

configs.load.id=[];             % file id numbers to use for analysis
configs.load.mincount=1000;
configs.load.maxcount=Inf;

configs.load.window{1}=[0.38,0.44];         % T [s]
configs.load.window{2}=[-35e-3,35e-3];      % X [m]
configs.load.window{3}=[-35e-3,35e-3];      % Y [m]