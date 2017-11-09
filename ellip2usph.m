function zxy_usph = ellip2usph(zxy,ecent,erad,evec,verbose)
% ZXY_USPH = ELLIP2USPH(ZXY, ECENT, ERAD, EVECS, VERBOSE)
%
% ellipsoid to unit-sphere 3D transform of vectors 
%
% see ellipsoid_fit.m for the ellipsoid parameters e___
% specifically, ellip params are defined in the XYZ coord sys
%
% ZXY is a Nx3 array of cartesian vectors (ARRAY)
%

nvecs=size(zxy,1);
xyz_temp=circshift(zxy,-1,2);

% centre to ellipsoid
% TODO - ecent seems to be a column array: a row vect would be better
xyz_temp=xyz_temp-repmat(ecent',nvecs,1);

% rotate to ellipsoid principal axis coords
M_rot=evec;   % rotation matrix: X'Y'Z'(principal ellipsoid) --> XYZ
xyz_temp=(M_rot\xyz_temp')';    % inverse transform

% map to unit-sphere by scaling along principal axis
% TODO - erad seems to be a column array: a row vect would be better
xyz_temp=xyz_temp./repmat(erad',nvecs,1);

% Rotate back to original (detector) axis
xyz_temp=(M_rot*xyz_temp')';

% we work in the ZXY coord system
zxy_usph=circshift(xyz_temp,1,2);

%%% TODO verbose output to show the before/after vecs

end