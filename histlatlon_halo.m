function [n] = histlatlon_halo(v_zxy,latlon_ed)
%   Histogram for halo data (Cart vecs in ZXY system) in the 2D lat-lon
%   grid.
%
%   latlon_ed: 1x2 cell-array of bin edges in lat/long-itude
%   

% Cart to sph-polar
vs=zxy2sphpol(v_zxy);

% 2D histogram
n=nhist(vs(:,1:2),latlon_ed);

end