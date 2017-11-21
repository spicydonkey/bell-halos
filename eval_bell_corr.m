%% Zonal analysis
%%% config zones, binning
nazim=configs.zone.nazim;
nelev=configs.zone.nelev;

binmethod=configs.zone.binmethod;
binwidth=configs.zone.binwidth;

%%% build grid of zones
% TODO - the closed sphere problem
switch binmethod
    case 1
        azim_vec=linspace(-pi,pi,nazim);  % edges

        % scan over all elev angle
        elev_vec=linspace(-pi/2,pi/2,nelev);    % easier to set nelev
        
        azim_cent=azim_vec;
        elev_cent=elev_vec;

    case 2
        error('binmethod=2 is unavailable.');
end
[azim_grid,elev_grid]=ndgrid(azim_cent,elev_cent);

%%% get counts in zone
% TODO - currently binning with fixed bin width - try Gaussian convolution
nn_halo=cell(2,1);
for ii=1:2
    nn_halo{ii}=halo_zone_density(halo_k0(:,ii),azim_vec,elev_vec,binwidth,binmethod);
    nn_halo{ii}=cat(3,nn_halo{ii}{:});  % nazim x nelev x nshot array
end

%%% statistics
nn_halo_mean=cellfun(@(x)mean(x,3),nn_halo,'UniformOutput',false);
nn_halo_std=cellfun(@(x)std(x,0,3),nn_halo,'UniformOutput',false);
nn_halo_serr=cellfun(@(x)std(x,0,3)/sqrt(size(x,3)),nn_halo,'UniformOutput',false);

% %% plot - zone count distribution
% hfig_halo_zone_density=figure();
% for ii=1:2
%     subplot(1,2,ii);
%     plot_sph_surf(azim_grid,elev_grid,nn_halo_mean{ii});
%     axis on;
%     box on;
%     xlabel('$K_X$');
%     ylabel('$K_Y$');
%     zlabel('$K_Z$');
%     title(configs.halo{ii}.string);
%     c=colorbar('SouthOutside');
%     c.Label.String='Avg. counts';
%     c.TickLabelInterpreter='latex';
%     c.Label.Interpreter='latex';
% end

%% cull background streak
% mf=1 halo looks fine
max_nn=max(max(nn_halo_mean{2}));   % max counts in zone from mf=1

% background from spontaneous scattering leaves a bright streak in mf=0
nn_thresh=max_nn*2;     % TODO - seems to work
bool_thresh=(nn_halo_mean{1}>nn_thresh);    % boolean on grid to treat as noisy

%%% cull mf=0
% set over threshold to NaN
nn_halo{1}(repmat(bool_thresh,[1,1,size(nn_halo{1},3)]))=NaN;

% redo statistics - shouldn't change!
%%% statistics
nn_halo_mean=cellfun(@(x)mean(x,3),nn_halo,'UniformOutput',false);
nn_halo_std=cellfun(@(x)std(x,0,3),nn_halo,'UniformOutput',false);
nn_halo_serr=cellfun(@(x)std(x,0,3)/sqrt(size(x,3)),nn_halo,'UniformOutput',false);

% %% plot - background filtered distribution
% hfig_halo_zone_density=figure();
% for ii=1:2
%     subplot(1,2,ii);
%     plot_sph_surf(azim_grid,elev_grid,nn_halo_mean{ii});
%     axis on;
%     box on;
%     xlabel('$K_X$');
%     ylabel('$K_Y$');
%     zlabel('$K_Z$');
%     title(configs.halo{ii}.string);
%     c=colorbar('SouthOutside');
%     c.Label.String='Avg. counts';
%     c.TickLabelInterpreter='latex';
%     c.Label.Interpreter='latex';
% end

%% Mixing characterisation
% Rabi cycle coefficient
P_rabi=nn_halo_mean{2}./(nn_halo_mean{1}+nn_halo_mean{2});
% 1 when all state in mf=1; 0 when all states in mf=0

% % plot
% hfig_rabi_coeff=figure();
% plot_sph_surf(azim_grid,elev_grid,P_rabi);
% 
% axis on;
% box on;
% xlabel('$K_X$');
% ylabel('$K_Y$');
% zlabel('$K_Z$');
% % title(configs.halo{ii}.string);
% c=colorbar('SouthOutside');
% c.TickLabelInterpreter='latex';
% c.Label.Interpreter='latex';
% % c.Label.String='$P(m_F=1)$';
% c.Label.String='$P(\uparrow)$';

%% Bell correlations
% TODO - fix all uncertainties in this!

nn_sum=nn_halo{1}+nn_halo{2};   % total number of atoms with same momentum
Jz=nn_halo{2}-nn_halo{1};       % total spin of atoms with same momentum

clearvars nn_halo;

% simple transform to match back-to-back momentum modes
nn_sum_bb=flip_bb(nn_sum,1);
Jz_bb=flip_bb(Jz,1);

% evaluate correlation coefficient
JJ=Jz.*Jz_bb;

clearvars Jz Jz_bb;

NN=nn_sum.*nn_sum_bb;

clearvars nn_sum nn_sum_bb;

E=mean(JJ,3)./mean(NN,3);

clearvars JJ NN;

% TODO
% each grid point has an associated (theta,phi) rotation angle
% Bell inequality in CHSH form is a linear combination of 4 correlation
% coefficients at a precisely related rotation angle configuration.

% %% plot correlations
% hfig_corr=figure();
% plot_sph_surf(azim_grid,elev_grid,E);
% 
% axis on;
% box on;
% xlabel('$K_X$');
% ylabel('$K_Y$');
% zlabel('$K_Z$');
% % title(configs.halo{ii}.string);
% c=colorbar('SouthOutside');
% c.TickLabelInterpreter='latex';
% c.Label.Interpreter='latex';
% c.Label.String='$E(\theta,\phi)$';
