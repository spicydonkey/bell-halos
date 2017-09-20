% Analyse correlations from Bell test
%
%

%% config
path_data_dir='C:\Users\HE BEC\bell';           % path to data directory
path_data_bell='bell_v1_3_20170920_20_15_tight';        % bell test
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

%%% rotation angles
% get index to Raman mixing amp
idx_Theta=cellfun(@(x) find(x==raman_amp),S_loop.ampRaman_mf);

% get rotation angles at different momenta (local operation)
Theta=zeros([size(E),2]);       % ELEV X AZIM X MF
for ii=1:2
    Theta(:,:,ii)=S_loop.Theta{ii}(:,:,idx_Theta(ii));
end

% NOTE - mf=0 Theta is wrong - use mf=1 result from symmetry
Theta(:,:,1)=Theta(:,:,2);      % assume mf=0 rotate by approx equal angle to mf=1 for equal momentum
% Theta(:,:,1)=-Theta(:,:,2);

% difference in rotation angle between back-to-back zones
dTh_bb=Theta(:,:,2)-flip_bb(Theta(:,:,1));

%% 
%%% plot correlation vs theta
hfig_corr_v_theta=figure(1);
plot(dTh_bb(:),E(:),'o');

xlim([0,pi]);
ylim([-1,1]);
xlabel('$\Delta \theta$');
ylabel('$E$');

% hidden variable prediction
hold on;
line([0,pi],[-1,1],'Color','r');