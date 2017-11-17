% Real-time monitor for mode occupancy
%
%

% USER CONFIGS
verbose=1;      % 0 / 1 

dir_data='\\amplpc29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output';

cfg.txy_window={[0.38,0.44],[-35,35]*1e-3,[-35,35]*1e-3};
cfg.bec_pos_1=[1.6524,-3.3e-3,4.6e-3];
cfg.bec_pos_2=[1.7055,-2.8e-3,4.8e-3];
cfg.bec_r=8e-3;
cfg.r_th=bec_r;
cfg.halo_dR=0.1;        % mf=0 captures very well
cfg.elev_max=asin(0.8);

nShotsLatest=50;

% END USER CONFIGS

while true
    close all;

    % run it
    this_mocc=getLatestModeOcc(dir_data,nShotsLatest,cfg,verbose);
    
    % store it in a ring buffer and plot it


    % waiting for some more data...
    pause('on')
    pause(30)      % it doensn't need to run every shot
end


