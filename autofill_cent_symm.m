function [azf,elf,ff]=autofill_cent_symm(az0,el0,f0)
% AUTOFILL_CENT_SYMM fills central symmetric distribution given 
% a hemisphere region az-elev (2D) surface distribution
% 
% * az0/el0, are 2D grid of azim/elev angles in radian.
% * DIM 1 must be AZIM (DIM 2 == ELEV)
% * 
% * el0 must be symmetric around 0 (equator)
%
% DKS
% 2018-10-31
% 

% cycle azimuthally to fill theta in [-pi,0) and pad theta=pi (theta=0)
azf=[az0-pi;az0;az0(1,:)+pi];
elf=[el0;el0;el0(1,:)];
ff=[fliplr(f0);f0;fliplr(f0(1,:))];     % flip elev

end