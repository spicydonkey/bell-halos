% Misc postprocessing of pair-src raw data to create pretty density plot
%
% DKS
% 20180606

% TODO
%   [ ] density image - not indiv counts
%


%% config
%%% data
f_data='out_pairsource_20180606.mat';
path_data='C:\Users\HE BEC\bell\2018april\exp5_bell_yrot\postproc';

%%% filter
% circ filter
filt_cyl_cent=[0,1.4e-3,4.2e-3];
filt_cyl_r=31e-3;

% hotspot
filt_hotspot_cent=[0,30.7e-3,-0.5e-3];
filt_hotspot_r=2.5e-3;

% z-window
filt_zwin_lim=[1.61,1.76];


%%% vis
[c1,c2]=palette(10);
id_col=1;

mark_area=5;
mark_type='o';
mark_col=c2(id_col,:);
mark_alph=1;        % DON'T adjust!
line_col='none';
line_alph=1;        % DON'T adjust!
line_width=0.5; 

%% load data
S=load(fullfile(path_data,f_data));

zxy_raw=S.zxy_src;      % get raw zxy data of pair-src (orig halo)


% % vis
% figure('Name','raw');
% plot_zxy(zxy_raw);
% axis equal;
% xlabel('x');
% ylabel('y');
% zlabel('z');


%% filter
% collate all shots
zxy_filt=vertcat(zxy_raw{:});       

% cyl filter
zxy_filt=zxy_filt-filt_cyl_cent;

r_cyl=vnorm(zxy_filt(:,2:3),2);
b_cyl_filt=r_cyl<filt_cyl_r;

zxy_filt=zxy_filt(b_cyl_filt,:);


% hotspot filter
dr_hotspot=vnorm(zxy_filt(:,2:3)-filt_hotspot_cent(2:3),2);
b_hotspot_filt=dr_hotspot>filt_hotspot_r;

zxy_filt=zxy_filt(b_hotspot_filt,:);


% z-window filter
b_zwin_filt=(zxy_filt(:,1)>filt_zwin_lim(1))&(zxy_filt(:,1)<filt_zwin_lim(2));
zxy_filt=zxy_filt(b_zwin_filt,:);



%% vis
figure('Name','src_filt');
% plot_zxy({zxy_filt},[],10);
scatter3(zxy_filt(:,2),zxy_filt(:,3),zxy_filt(:,1),...
    mark_area,...
    'Marker',mark_type,...
    'MarkerFaceColor',mark_col,...
    'MarkerFaceAlpha',mark_alph,...
    'MarkerEdgeColor',line_col,...
    'MarkerEdgeAlpha',line_alph,...
    'LineWidth',line_width);

axis equal;
% axis off;
box on;

set(gca,'XTick',[]);
set(gca,'YTick',[]);
set(gca,'ZTick',[]);

xlabel('$x$'); 
ylabel('$y$');
zlabel('$z$');

set(gca,'FontSize',15);

