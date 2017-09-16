function halo_zone_density(halo,verbose)
% zonal analysis on scattered halo
% halo_zone_density(halo, verbose)
% 
% halo: num_shot x 1 cell-array of num_countsx3 double array (XYZ)
%

%% evaluate normalised density per unit solid angle at all angular zones
    
% XYZ --> spherical angle (2D)
% set up grid of spherical angles + bin style + width
    % get effective counts at each grid point
        % per shot
    % evaluate normalised density (integrates to 1 around domain)
% statistics
    % mean
    % sdev, serr (shot-to-shot)
    % bootstrapping
    
% return
    % grid
    % normalised density statistics at grid points