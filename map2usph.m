function [q,v_ellip3]=map2usph(k)
% simple transformation to ~unit sphere by the ellipsoid_fit routine
% q=map2usph(k)
% 
%   k: cell-array of ZXY vector-lists
%
%   v_ellip3: ellipsoid parameters incl center, radii, principal axes,
%       algebraic form params
%
% DKS
%

% collate dataset as a Nx3 array
K=vertcat(k{:});

% fit ellipsoid (XYZ coords)
[ec,er,ev,v]=ellipsoid_fit(zxy2xyz(K),'');     % in XYZ frame

% get fitted ellipsoid parameters
v_ellip3.cent=ec;   % center
v_ellip3.rad=er;    % radii
v_ellip3.vaxis=ev;  % principal axes
v_ellip3.v=v;       % algebraic param vector

% transform data to unit-sphere
q=cellfun(@(x)ellip2usph(x,ec,er,ev),k,'UniformOutput',false);

end
