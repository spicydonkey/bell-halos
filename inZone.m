function bool_inzone=inZone(th_1,phi_1,th_2,ph_2,dpsi)
% BOOL_INZONE = INZONE(AZIM,ELEV,THETA,PHI,DPSI)
%
% gets indices to spherical zones (defined by azim-elev grid) in the
% vicinity of reference direction
%

% get difference angle at all zones
psi=diffAngleSph(th_1,phi_1,th_2,ph_2);

% evaluate logical condition
bool_inzone=psi<dpsi;

end
