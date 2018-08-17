%% Fig2: Generate 2D-g2 profile of pair-source 
% DKS
% 2018-08-17
%% load data
path_data='C:\Users\HE BEC\Documents\lab\bell_momentumspin\pairsource\90\rev\run_4_rev.mat';
load(path_data);
%% plot and graphics
[h,s]=surf_g2_2d(g2{1}{3},dk{1},idx_dk0);   % surf plot
colormap magma;
s.EdgeColor='none';     % turn off edge
box off;
grid off;
%%% mesh
submeshsurf(s,0.5,'w',0.5,true);        % dk grid is an nd-type grid
%%% annotation
zlabel('$g^{(2)}_{\mathrm{BB},\uparrow\downarrow}(\Delta \mathbf{k})$');
ax=gca;
ax.FontSize=24;
h.Renderer='Painters';