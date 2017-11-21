% Generate mixing angle map for ideal local operation
%
%

% Set configs
if ~exist('OVERRIDE_CONFIG_FLAG','var')    
    % fun_localoper=@(th,phi) 0*phi;
    fun_localoper=@(th,phi) pi*sin(20*th).*sin(20*phi);
    
    verbose=1;
    
    nazim=100;
    nelev=50;
end

%%  
if ~exist('fun_localoper','var')
    error('fun_localoper undefined!');
end


% build azim / elev
avec=linspace(-pi,pi,nazim);            
evec=linspace(-pi/2,pi/2,nelev);        % pole-to-pole

[azim{1},elev{1}]=ndgrid(avec,evec);      % data typing conforms with main analysis

% proxy for simulated ideal mixing map
ampRaman_mf={NaN};      % NaN is for simulation


%% Evaluate rotation angle from given function
Theta=cell(2,1);
Theta{1}=fun_localoper(azim{1},elev{1});    % evaluate state rotation angle
Theta{2}=Theta{1};      % mf=1 rotates same as mf=0


%% plot map of rotation angle
if verbose>0
    hfig_ideal_theta=figure();
    for ii=1:2
        subplot(1,2,ii);
        plot_sph_surf(azim{1},elev{1},Theta{ii});
        
        axis on;
        xlabel('$K_x$');
        ylabel('$K_y$');
        zlabel('$K_z$');
        cbar=colorbar('SouthOutside');
        cbar.Label.String='$\theta$';
        
        title(sprintf('$m_F=$%d',ii-1))
        cbar.TickLabelInterpreter='latex';
        cbar.Label.Interpreter='latex';
    end
end