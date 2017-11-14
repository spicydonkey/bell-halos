function Fbb_thphi=flip_bb_new(F_thphi,dim_azim)
%
% mat_thphi: 2D matrix of values defined on wrapped spherical grid
%   elev x azim mesh (Note: from meshgrid(azim,elev))
%   elev angles should be symmetric; azim completely wraps 2*pi
% mat_flipped: angle inverted (back-to-back on sphere)

% Invert spherical map about the origin
%
% FBB_THPHI = FLIP_BB_NEW(F_THPHI,DIM_AZIM)
%
% F_THPHI: 2D matrix of values defined on a spherical grid given as:
%   - grid: azim x elev
%   - elev angles MUST be symmetric
%   - azim angles MUST be complete `modulo 2pi` but UNIQUE!
%
% DIM_AZIM: `azimuthal` dimension
%   - ndgrid / ngrid    --> 1
%   - meshgrid / psycho --> 2 (WHY WOULD YOU EVEN USE THIS?! HUH?)
%
% TODO
%   - is there ideal dimension length parity?
%
% DKS 
% 2017/11/14
%

% get the complementary azimuthal dimension
dim_elev=mod(dim_azim,2)+1;

% ROTATE the azimuthal by 180 deg; FLIP (mirror) the polar (elev)
nazim=size(F_thphi,dim_azim);
if mod(nazim,2)~=0
    warning('Azimuthal length is NOT EVEN. Spherical inversion will not be exact.');
end
n_rot_azim=floor(nazim/2);
% should REQUIRE EVEN azim length array - this is the ONLY way to for 
% filp_bb to be involutory (self-inverse)!

Fbb_thphi=circshift(F_thphi,n_rot_azim,dim_azim);  % rotate azimuthal angle by approx pi
Fbb_thphi=flip(Fbb_thphi,dim_elev);        % flip elevation

end