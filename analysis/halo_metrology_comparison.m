% common
qe = 0.1;


%% Magnetometry
magnetometry.sig_B = 4e-3;       % magnetic field uncertainty (G)
magnetometry.dx = 360e-6;           % halo diameter
magnetometry.tau_shot_max = 3.55e-6;    % max interrogation time per shot (s)
magnetometry.n_tau = 32;                % # interrogation times
magnetometry.shot_per_tau = 10;         % # shots per tau
magnetometry.bin_alpha = 0.0619*pi;     % bin half-cone angle (rad)
magnetometry.n_in_bin = 64;             % avg. # atoms in bin


magnetometry.sig_dBdx = sqrt(2)*magnetometry.sig_B/magnetometry.dx;
magnetometry.solid_angle = cone_solid_angle(magnetometry.bin_alpha);
magnetometry.N_tot_bin = magnetometry.n_in_bin * magnetometry.n_tau * magnetometry.shot_per_tau;

magnetometry.sig_dBdx*(sqrt(magnetometry.N_tot_bin)...
    *magnetometry.tau_shot_max)


%% Gradiometry
gradiometry.sig_dBdx = 0.2;   % G/m
gradiometry.tau_shot_max = 1.7e-3;
gradiometry.bin_alpha = 0.124*pi;
gradiometry.N_tot_bin = 1770;


gradiometry.solid_angle = 2*cone_solid_angle(gradiometry.bin_alpha);

gradiometry.sig_dBdx*(sqrt(gradiometry.N_tot_bin)...
    *gradiometry.tau_shot_max)


%% functions
function solid_angle = cone_solid_angle(alpha)
    solid_angle = 2*pi*(1 - cos(alpha));
end