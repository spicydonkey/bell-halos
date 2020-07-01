function E = mdl_rhop_offres(beta,x)
%mdl_rhop_offres is a E(th,th) model for rho+ state (realistic Bell triplet)
%with off-resonant rotation.
%
%   BETA is a 2-length vector: [p, delta]
%       p       :   mixing parameter, 0 <= p <= 1
%       (d)elta :   detuning parameter (elev angle of rotation axis):
%       |delta| <= pi/2
%
%   The model form is:
%   y = -p * (sin(d)^2 + cos(d)^2*cos(x)).^2 ...
%       + p * ( (sin(d)*cos(d)*(1-cos(x))).^2 + (cos(d)*sin(x)).^2)
%
%   DKS 2020

p = beta(1);            
delta = beta(2);

E = -p * (sin(delta)^2 + cos(delta)^2*cos(x)).^2 ...
      + p * ( (sin(delta)*cos(delta)*(1-cos(x))).^2 + (cos(delta)*sin(x)).^2);
end