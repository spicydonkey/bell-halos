function psi=sphdiffangle(th1,phi1,th2,phi2)
% psi = sphdiffangle(th1, phi1, th2, phi2)
%
% th: azimuthal angles measured from X-axis in the positive Z-rot direction
% phi: elev angles measured from XY-plane
% 

psi=acos(sin(phi1).*sin(phi2)+cos(phi1).*cos(phi2).*cos(th1-th2));
end