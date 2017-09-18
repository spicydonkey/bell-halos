function nn=halo_zone_density(halo,azim,elev,binwidth,verbose)
% zonal analysis on scattered halo
% halo_zone_density(halo, verbose)
% 
% halo: num_shot x 1 cell-array of num_countsx3 double array (ZXY)
    % halo should be correctly centered
% azim: def grid points - array of azim angles
% elev: def grid points - array of elev angles (meas from XY plane)
%
% nn: num_shot x 1 cello-array of meshgrid of double-array of normalised number density

%% evaluate normalised density per unit solid angle at all angular zones
% TODO
    % meshgrid for azim-elev is ideal for surf plotting on sphere

% ZXY --> spherical angle (2D)
nshot=size(halo,1);
halo_sphpol=cell(nshot,1);
for ii=1:nshot
    this_halo=halo{ii};
    halo_sphpol{ii}=zeros(size(this_halo,1),2);
    % note: halo is defined in ZXY coord system
    [halo_sphpol{ii}(:,1),halo_sphpol{ii}(:,2)]=cart2sph(this_halo(:,2),this_halo(:,3),this_halo(:,1));
end

grid_thphi=[azim(:),elev(:)];
ngrid=size(grid_thphi,1);       % number of grid points
N_zone=cell(nshot,1);           % counts in zones per shot

% method 1
% bin style: simple counting bin with limits
% binwidth=deg2rad(10);   % some appropriate binning angle (rad)
for ii=1:nshot
    % set up grid of spherical angles + bin style + width
    % get effective counts at each grid point
    % per shot
    this_N_zone=zeros(size(azim));   % preallocate counts in zones
    
    this_halo_thphi=halo_sphpol{ii};
%     this_ncounts=size(this_halo_thphi,1);
    
    % for each grid-point evaluate angle difference
    for jj=1:ngrid
        this_thphi=grid_thphi(jj,:);  % this grid point
        % evaluate diff angles to the grid point
        this_psi=sphdiffangle(this_thphi(1),this_thphi(2),this_halo_thphi(:,1),this_halo_thphi(:,2));
        
        % count number in this zone
        this_N_zone(jj)=sum(this_psi<binwidth);
    end
    N_zone{ii}=this_N_zone;
end

% evaluate normalised density (integrates to 1 around domain)
    % NOTE: simple linspace sampling for azim,elev from limits doesn't
    % uniformly sample the sphere.

% statistics
    % mean
    % sdev, serr (shot-to-shot)
    % bootstrapping
    
% return
    % normalised density statistics at grid points
    
% return the shot-wise counts in zone
nn=N_zone;