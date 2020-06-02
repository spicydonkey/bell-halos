function [n_r,r_cent]=rdist_truncsph(r,r_ed,z_trunc)
%radial distribution in truncated sphere
%
%   r:          array of vector norms
%   r_ed:       radial bin edges
%   z_trunc:    truncation z from equator
%
%   n_r: radial density profile (per volume)
%   r_cent: bin centers
%

% histogram counts
N_r=histcounts(r,r_ed);

% get bin centers
r_cent=edge2cent(r_ed);

% normalise (TRUE to preproc data) by truncated spherical shell
vol_r_bin = vol_truncsphshell(r_ed(1:end-1),r_ed(2:end),z_trunc);

n_r=N_r./vol_r_bin;                 % normalise by bin vol
end