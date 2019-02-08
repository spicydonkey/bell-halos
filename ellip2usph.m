function zxy_usph = ellip2usph(zxy,ecent,erad,evec)
% ZXY_USPH = ELLIP2USPH(ZXY, ECENT, ERAD, EVECS, VERBOSE)
%
% ellipsoid to unit-sphere 3D transform of vectors 
%
% see ellipsoid_fit.m for the parameterisation of fitted elllipsoid
% specifically, ellip params are defined in the XYZ coord sys
%
% ZXY is a Nx3 array of cartesian vectors (ARRAY)
%

% to XYZ coord sys
xyz_temp=zxy2xyz(zxy);

% centre to ellipsoid
xyz_temp=xyz_temp-ecent';

% rotate to ellipsoid principal axis coords
M_rot=evec;         % rotation matrix: X'Y'Z'(principal ellipsoid) --> XYZ
xyz_temp=(M_rot\xyz_temp')';    % inverse transform (and keep vectors as rows)

% map to unit-sphere by scaling along principal axis
xyz_temp=xyz_temp./(erad');

% Rotate back to original (detector) axis
xyz_temp=(M_rot*xyz_temp')';

% to ZXY coord system
zxy_usph=xyz2zxy(xyz_temp);

end