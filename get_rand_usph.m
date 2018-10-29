function v = get_rand_usph(N)    
% Generate a randomly sampled point on a 3D unit sphere
%   V = get_rand_usph(N)
%   
%   N:  number of points
%   v:  Nx3 array of cartesian 3-vectors
%
% DKS
% 2018-10-29
%
    theta=2*pi*rand(N,1);           % azim (rad)
    sin_phi=2*rand(N,1)-1;          % sin(phi)
    cos_phi=sqrt(1-sin_phi.^2);     % cos(phi)
    
    v=[cos_phi.*cos(theta),cos_phi.*sin(theta),sin_phi];
end