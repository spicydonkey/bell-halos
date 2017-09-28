% Generate a randomly sampled point on a 3D unit sphere
function u_xyz = get_rand_usph(N)
    rand_th_azim=rand(N,1);
    rand_th_elev=rand(N,1);
    
    th_azim=2*pi*rand_th_azim;
    th_elev=asin(2*rand_th_elev-1);
    
    u_xyz=[cos(th_elev).*cos(th_azim),cos(th_elev).*sin(th_azim),sin(th_elev)];
end