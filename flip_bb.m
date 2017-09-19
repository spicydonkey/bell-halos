function mat_flipped=flip_bb(mat_thphi)
% mat_flipped = flip_bb(mat_thphi)
%
% mat_thphi: 2D matrix of values defined on wrapped spherical grid
%   elev x azim mesh (Note: from meshgrid(azim,elev))
%   elev angles should be symmetric; azim completely wraps 2*pi
% mat_flipped: angle inverted (back-to-back on sphere)

nazim=size(mat_thphi,2);
n_rot_azim=floor(nazim/2);

mat_flipped=circshift(mat_thphi,n_rot_azim,2);  % rotate azimuthal angle by approx pi
mat_flipped=flipud(mat_flipped);        % flip elevation

end