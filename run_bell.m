% Analyse correlations from Bell test
%
%

%% config
path_data_dir='C:\Users\HE BEC\bell';           % path to data directory
path_data_bell='idealsource_20171003_loop_test3';        % bell test
path_data_loop='idealloop_20171003_test3';          % local operations

% raman_amp=0.37;     % this run's LO config (raman amplitude)
% raman_amp=0;
raman_amp=NaN;      % for simulated local mixing


%% load data
% load Bell test results (incl. experimental correlations)
S_bell=load(fullfile(path_data_dir,path_data_bell));

% load local operations
S_loop=load(fullfile(path_data_dir,path_data_loop));

%% 
% get zone
azim_g=S_bell.azim_grid;
elev_g=S_bell.elev_grid;

% get correlations
E=S_bell.E;

%% rotation angles
% get index to Raman mixing amp
if isnan(raman_amp)
    idx_Theta=[1,1];
else
    idx_Theta=cellfun(@(x) find(x==raman_amp),S_loop.ampRaman_mf);
end

% get raw rotation angles at different momenta (local operation)
azim_raw=S_loop.azim{1};    % TODO - we know azim{ii} are the same
elev_raw=S_loop.elev{1};

Theta_raw=zeros([size(azim_raw),2]);       % ELEV X AZIM X MF
Theta_raw=[];      % rotation angle loaded from file (coarsely meshed?)
for ii=1:2
    Theta_raw(:,:,ii)=S_loop.Theta{ii}(:,:,idx_Theta(ii));
end

% NOTE - mf=0 Theta is wrong - use mf=1 result from symmetry
Theta_raw(:,:,1)=Theta_raw(:,:,2);      % assume mf=0 rotate by approx equal angle to mf=1 for equal momentum
% TODO - for debugging
if raman_amp==0
    % assume all rotation is zero - which it should be
    Theta_raw=zeros(size(Theta_raw));
end

% difference in rotation angle between back-to-back zones
dTh_bb_raw=Theta_raw(:,:,2)-flip_bb(Theta_raw(:,:,1));

% interpolate to build diff rotation angle at Bell data
% TODO - it's only approximate now - pretty messy!
if isnan(raman_amp)
    dTh_bb=interp2(azim_raw,elev_raw,dTh_bb_raw,azim_g(1:end-1,1:end-1),elev_g(1:end-1,1:end-1));
else
    dTh_bb=interp2(azim_raw(1:end-1,1:end-1),elev_raw(1:end-1,1:end-1),dTh_bb_raw,azim_g(1:end-1,1:end-1),elev_g(1:end-1,1:end-1));
end

% plot interpolated rotation angle
hfig_dTh_interp=figure();
plot_sph_surf(azim_g,elev_g,dTh_bb);
title('Interpolated rotation angles');
axis on;
xlabel('$K_x$');
ylabel('$K_y$');
zlabel('$K_z$');
cbar=colorbar('SouthOutside');
cbar.Label.String='$\Delta \theta$';

cbar.TickLabelInterpreter='latex';
cbar.Label.Interpreter='latex';

%% Process data
%%% cull elev edges
% NOTE: culling is as simple as setting those pixel values by NaN
n_elev_edge_cull=1;     % cull top and bottom by this many pixels
E(1:n_elev_edge_cull,:)=NaN;
E(end-n_elev_edge_cull+1:end,:)=NaN;

%%% Sort
% sort data
EE=E(:);
DTh=dTh_bb(:);

% remove NaNs
DTh=DTh(~isnan(EE));
EE=EE(~isnan(EE));
EE=EE(~isnan(DTh));
DTh=DTh(~isnan(DTh));

[DTh,Isort]=sort(DTh,'ascend');
EE=EE(Isort);

%%% RAW plot: correlation vs theta
hfig_corr_v_theta_raw=figure();
plot(DTh,EE,'bo','MarkerSize',3);

% xlim([0,pi]);
xlim([-pi,pi]);
% ylim([-1,1]);
xlabel('$\Delta \theta$');
ylabel('$E$');

% hidden variable prediction
hold on;
line([-pi,0],[1,-1],'Color','r','LineWidth',1.5);
line([0,pi],[-1,1],'Color','r','LineWidth',1.5);
hold off;
box on;

%% smooth correlations
% Dth is "averaged" angular difference (vs. DTh)

n_Dth_bin=21;
Dth_edge=linspace(-pi,pi,n_Dth_bin+1);
Dth_cent=Dth_edge(1:end-1)+0.5*diff(Dth_edge);

Dth=zeros(2,n_Dth_bin);   % mean and serr
Eth=zeros(2,n_Dth_bin);  

for ii=1:n_Dth_bin
    bool_bin=(DTh>Dth_edge(ii)&DTh<Dth_edge(ii+1));     % boolean for items in bin
    n_items=sum(bool_bin);      % number of items in this bin
    
    % get values
    DTh_bin=DTh(bool_bin);
    EE_bin=EE(bool_bin);
    
    % statistics
    Dth(1,ii)=mean(DTh_bin);
    Dth(2,ii)=std(DTh_bin)/sqrt(n_items);
    
    Eth(1,ii)=mean(EE_bin);
    Eth(2,ii)=std(EE_bin)/sqrt(n_items);
end

% DEBUG
disp(Eth(~isnan(Eth)));

%%% SMOOTH plot: correlation vs theta
% misc
markersize=5;
linewidth=1.5;
gray_col=0.5*ones(1,3);         % gray data points
namearray={'LineWidth','MarkerFaceColor','Color'};      % error bar graphics properties
valarray={linewidth,'w','b'};                 % 90 deg (normal) data

hfig_corr_v_theta_smooth=figure();

%%% DATA
hcorr=ploterr(Dth(1,:),Eth(1,:),Dth(2,:),Eth(2,:),'o','hhxy',0);
set(hcorr(1),namearray,valarray,'MarkerSize',markersize,'DisplayName','Experiment');	% DATAPOINT
set(hcorr(2),namearray,valarray,'DisplayName','');                  % Y-err
set(hcorr(3),namearray,valarray,'DisplayName','');                  % X-err

%%% THEORIES
% LHV
hold on;
hLHV=line([-pi,0,pi],[1,-1,1],'LineStyle','--','Color',gray_col,'LineWidth',linewidth,'DisplayName','LHV');
hold off;

% QM
hold on;
dth_qm=linspace(-pi,pi,100);
E_qm=-cos(dth_qm);
hQM=plot(dth_qm,E_qm,'LineStyle','-','Color','k','LineWidth',linewidth,'DisplayName','QM');
hold off;

% annotate
ax=gca;

uistack(hQM,'bottom');
uistack(hLHV,'bottom');
oleg=legend([hcorr(1),hQM,hLHV],'Location','southeast');
xlabel('$\Delta \theta = \theta_1 - \theta_0$');
ylabel('$E$');
xlim([-pi,pi]);
ylim([-1,1]);

oleg.FontSize=11;
ax.FontSize=13;