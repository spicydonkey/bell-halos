function bool_inzone=inZone(azim,elev,theta,phi,dpsi)
% BOOL_INZONE = INZONE(AZIM,ELEV,THETA,PHI,DPSI)
%
% gets indices to spherical zones (defined by azim-elev grid) in the
% vicinity of reference direction
%

% get difference angle at all zones
psi=diffAngleSph(azim,elev,theta,phi);

% evaluate logical condition
bool_inzone=psi<dpsi;

end
