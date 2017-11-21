% Real-time monitor for mode occupancy
%
%

%-------------------------------------------------
%%% USER CONFIGS
verbose=0;      % 0 (just the trend) / 1 (message to stdout) / 2<= (all hell breaks loose. TODO silent figure)

dir_data='\\amplpc29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output';

cfg.txy_window={[0.38,0.44],[-35,35]*1e-3,[-35,35]*1e-3};
cfg.bec_pos_1=[1.6524,-3.3e-3,4.6e-3];
cfg.bec_pos_2=[1.7055,-2.8e-3,4.8e-3];
cfg.bec_r=8e-3;
cfg.r_th=cfg.bec_r;
cfg.halo_dR=0.1;        % mf=0 captures very well
cfg.elev_max=asin(0.8);

nShotsLatest=50;

t_pause=25;     % pretty good
lenCircBuff=500;

%%% END USER CONFIGS
%-------------------------------------------------

% figure with dynamic data source
hfig_trend=figure();
mocc_latest=NaN(lenCircBuff,1);
h=plot(mocc_latest,'YDataSource','mocc_latest',...
    'Color','k','LineStyle','--','Marker','o','LineWidth',2);
title(sprintf('Real-time mode occupancy monitor: %s',datetime));
ylabel('Mode occupancy');

while true
    %%%% tidy up
    clc;
    %figures to keep
    figs2keep = [hfig_trend];    
    all_figs = findobj(0, 'type', 'figure');
    delete(setdiff(all_figs, figs2keep));

    % run it
    this_mocc=getLatestModeOcc(dir_data,nShotsLatest,cfg,verbose);
    
    % store it in a ring buffer and plot it
    mocc_latest=circshift(mocc_latest,-1,1);
    mocc_latest(end)=this_mocc;

    % refresh plot
	title(sprintf('Real-time mode occupancy monitor: %s',datetime));
    refreshdata(h, 'caller')
    axis auto;
    drawnow;

    % waiting for some more data...
    pause('on')
    pause(t_pause)      % it doensn't need to run every shot
end


