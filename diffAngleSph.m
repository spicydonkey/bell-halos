function psi=diffAngleSph(th1,phi1,th2,phi2)
% PSI = SPHDIFFANGLE(TH1,PHI1,TH2,PHI2)
% Evaluates the relative angle between two vectors defined in spherical
% coordinates
% 
%
% TH: azimuthal angles (std Matlab convention)
% PHI: elev angles
% 

psi=acos(sin(phi1).*sin(phi2)+cos(phi1).*cos(phi2).*cos(th1-th2));
end