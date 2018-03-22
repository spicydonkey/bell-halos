function zxy_out = rad_deform(zxy,az,el,rr)
%Radially deform 3D vectors by a scaling map defined on the 2-sphere
%
%   zxy: cart vecs
%   az: azim grid
%   el: elev grid
%   rr: radial scaling factor map
%
%   zxy_out: Cart 3D vecs, scaled
%       vecs in ill-defined regions in rr are culled
%

% cart --> sph-pol
spol=zxy2sphpol(zxy);

% TODO - can do cleaner by a continuous map by fitting in sph-harms
% 0-2pi close discrete 2-sphere map
% NOTE: az-el should be standard sph-grid s.t. azims in DIM-1, elev in DIM-2.
az(end+1,:)=pi;
el(end+1,:)=el(end,:);
rr(end+1,:)=rr(1,:);

% scale norms
% NOTE: if map has NaNs then counts in those regions are culled
c=interpn(az,el,rr,spol(:,1),spol(:,2));    % interpolated radial scale factors
spol(:,3)=spol(:,3).*c;

% cull counts in ill-defined regions in given map
b_oor=isnan(c);
if sum(b_oor)>0
    spol=spol(~b_oor,:);
end

% sph-pol --> cart
zxy_out=sphpol2zxy(spol);

end