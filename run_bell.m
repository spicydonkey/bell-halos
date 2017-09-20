% Analyse correlations from Bell test
%
%

%% config
path_data_dir='C:\Users\HE BEC\bell';           % path to data directory
path_data_bell='bell_v1_3_20170920_wbb';        % bell test
path_data_loop='loop_v1_20170920_20_15';          % local operations

raman_amp=0.37;     % this run's LO config (raman amplitude)

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
idx_Theta=cellfun(@(x) find(x==raman_amp),S_loop.ampRaman_mf);

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

% difference in rotation angle between back-to-back zones
dTh_bb_raw=Theta_raw(:,:,2)-flip_bb(Theta_raw(:,:,1));

% interpolate to build diff rotation angle at Bell data
% TODO - it's only approximate now - pretty messy!
dTh_bb=interp2(azim_raw(1:end-1,1:end-1),elev_raw(1:end-1,1:end-1),dTh_bb_raw,azim_g(1:end-1,1:end-1),elev_g(1:end-1,1:end-1));

% plot interpolated rotation angle
hfig_dTh_interp=figure();
plot_sph_surf(azim_g,elev_g,dTh_bb);


%% process data
% sort data
EE=E(:);
DTh=dTh_bb(:);

% remove NaNs
DTh=DTh(~isnan(EE));
EE=EE(~isnan(EE));

[DTh,Isort]=sort(DTh,'ascend');
EE=EE(Isort);

%%% smoothing
% EE=smooth(EE,30);
% % gaussian convolution
% sigma=5;
% sz=10;
% x = linspace(-sz/2,sz/2,sz);
% gaussFilter=exp(-x.^2/(2*sigma^2));
% gaussFilter=gaussFilter/sum(gaussFilter); % normalize
% EE = conv(EE,gaussFilter,'same');

%%% plot correlation vs theta
hfig_corr_v_theta=figure(1);
plot(DTh,EE,'o');

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