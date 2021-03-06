function nn=halo_zone_density(halo,azim_vec,elev_vec,binwidth,binmethod,verbose)
% zonal analysis on scattered halo
% halo_zone_density(halo, verbose)
% 
% halo: num_shot x 1 cell-array of num_countsx3 double array (ZXY)
    % halo should be correctly centered
% azim_vec: def grid points - array of azim angles
% elev_vec: def grid points - array of elev angles (meas from XY plane)
%
% binwidth: width of bin
% binmethod:
%   1: simple - divide sphere into simple nonuniform distributed segments
%   (one-to-one). Note: binwidth will not unused.
%   2: angular - require binwidth (may be one-to-many OR miss counts).
%   
%
% nn: num_shot x 1 cell-array of normalised number density defined on spherical grid

%% parse inputs
if ~exist('binwidth','var')||~exist('binmethod','var')
    binwidth=[];        % Default: empty binwidth and set binmethod to simple
    binmethod=1;
end
if ~exist('verbose','var')
    verbose=0;
end

%% evaluate normalised density per unit solid angle at all angular zones
% ZXY --> spherical angle (2D)
nshot=size(halo,1);
halo_sphpol=cell(nshot,1);
for ii=1:nshot
    this_halo=halo{ii};
    halo_sphpol{ii}=zeros(size(this_halo,1),2);
    % note: halo is defined in ZXY coord system
    [halo_sphpol{ii}(:,1),halo_sphpol{ii}(:,2)]=cart2sph(this_halo(:,2),this_halo(:,3),this_halo(:,1));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % azim in range (-pi,pi]
    % elev in range [-pi/2,pi/2]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% build spherical grid zones
[azim_grid,elev_grid]=ndgrid(azim_vec,elev_vec);
nazim=length(azim_vec);
nelev=length(elev_vec);

% build edges
thphi_edge={azim_vec,elev_vec};     % see EDGE argument in "nhist.m"

grid_thphi=[azim_grid(:),elev_grid(:)];
ngrid=size(grid_thphi,1);       % number of grid points

N_zone=cell(nshot,1);           % counts in zones per shot

% evaluate count/intensity in each shot at each zone
for ii=1:nshot
    % set up grid of spherical angles + bin style + width
    % get effective counts at each grid point per shot
    
    this_halo_thphi=halo_sphpol{ii};
    
    % histogramming methods
    switch binmethod
        % build this_N_zone - counts/intensity at the defined spherical grid
        case 1
            % Simple spherical Lat-Lon grid 
            this_N_zone=nhist(this_halo_thphi,thphi_edge);     % 2D histogram - transpose
            
        case 2
            % Difference angle
            % NOTE: this method isn't one-to-one!
            % for each grid-point evaluate angle difference
            this_N_zone=zeros(nazim,nelev);   % preallocate counts in zones
            for jj=1:ngrid
                this_thphi=grid_thphi(jj,:);  % this grid point
                % evaluate diff angles to the grid point
                this_psi=sphdiffangle(this_thphi(1),this_thphi(2),this_halo_thphi(:,1),this_halo_thphi(:,2));
                
                % count number in this zone
                this_N_zone(jj)=sum(this_psi<binwidth);
            end
    end
    N_zone{ii}=this_N_zone;
end

% evaluate normalised density (integrates to 1 around domain)
    % NOTE: simple linspace sampling for azim,elev from limits doesn't
    % uniformly sample the sphere.

%%% statistics
% easier done outside
    % mean
    % sdev, serr (shot-to-shot)
    % bootstrapping
    
% return
    %? normalised density statistics at grid points
    
% return the shot-wise counts in zone
nn=N_zone;
