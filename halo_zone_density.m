function nn=halo_zone_density(halo,azim,elev,verbose)
% zonal analysis on scattered halo
% halo_zone_density(halo, verbose)
% 
% halo: num_shot x 1 cell-array of num_countsx3 double array (XYZ)
% azim: def grid points - array of azim angles
% elev: def grid points - array of elev angles (meas from XY plane)
%
% nn: meshgrid of double-array of normalised number density

%% evaluate normalised density per unit solid angle at all angular zones
% TODO
    % meshgrid for azim-elev is ideal for surf plotting on sphere

% XYZ --> spherical angle (2D)
% set up grid of spherical angles + bin style + width
    % get effective counts at each grid point
        % per shot
    % evaluate normalised density (integrates to 1 around domain)
        % NOTE: simple linspace sampling for azim,elev from limits doesn't
        % uniformly sample the sphere.
% statistics
    % mean
    % sdev, serr (shot-to-shot)
    % bootstrapping
    
% return
    % grid
    % normalised density statistics at grid points