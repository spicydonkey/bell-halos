function [bool_inzone, n] =inZone(th_1,phi_1,th_2,ph_2,dpsi)
% [BOOL_INZONE, N] = INZONE(AZIM,ELEV,THETA,PHI,DPSI)
%
%   gets indices to spherical zones (defined by azim-elev grid) in the
%   vicinity of reference direction
%
%   n: returns number of atoms in zone
%


% get difference angle at all zones
psi=diffAngleSph(th_1,phi_1,th_2,ph_2);

% evaluate logical condition
bool_inzone=psi<dpsi;

end
