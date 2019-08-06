function [Y_lm,PHI,THETA] = sphharmY(l,m,phi,theta)
% sphharmY: calculates the spherical harmonics in real form
%
%   phi: azimuthal angle (from +x-axis CC around +z-axis
%   theta: polar angle (measured from +z-axis)
%
% see: https://en.wikipedia.org/wiki/Spherical_harmonics
%   
% DKS 2019

[PHI,THETA] = ndgrid(phi,theta);

% legendre polynomial in polar angle
P_lm = legendre(cos(theta));

if m > 0
    Y_lm = NaN;
elseif m < 0
    Y_lm = NaN;
else
    % m == 0
    Y_lm = sqrt((2*l + 1)/(4*pi))*P_lm(1,:).*ones(size(phi));
end

end