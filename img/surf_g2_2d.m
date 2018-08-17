function [h,s]=surf_g2_2d(g2,dk,idx_dk0)
%   g2 should be a 3-dim array defined on
%   dk a 1x3 cell array of Z,X,Y grid points
%   idx_dk0 is index of dim @ dk=0 in that component
%
%   Cart dims are in EXP order: Z,X,Y
%
%   NOTE: currently ZX-plane through Y at 0
%


h=figure('Name','g2_2d');

s=surf(dk{1}(:,:,idx_dk0),dk{2}(:,:,idx_dk0),g2(:,:,idx_dk0));

% shading interp;

axis square;

xlabel('$\Delta k _z$');
ylabel('$\Delta k _x$');
zlabel('$g^{(2)}$');

end