% Testing Bell correlation distribution in SOURCE
%

tstart=tic;

%% load data
load('src1_data_digest.mat');
% loads entangled source halos (halo_k0)
K=halo_k0;


%% configure
% scattered modes
nAz=100;        % num of azimuthal divisions in [-pi,pi)
nEl=50;         % num of elevation divisions in [-pi/2,pi/2]
% NOTE: nAz should be EVEN to ensure flip_bb.m correctly matches
% spherically opposite modes

% counting
% count_mode='gauss';
% sig_mode=[0.018,Inf];
% lim_mode=[3,Inf];

count_mode='simple';
sig_mode=0.1;
lim_mode=[];

z_max_nan=0.7;      % max-z used to filter halos as bad

% construct sph-grid zones
[Az,El]=sphgrid(nAz,nEl);


%% get counts in momenta modes
% TODO 
%   [ ] doesn't need to be run (?) for new halo k-space dataset if we could
%   manipulate this
%   [x] Parallelise
%       1. serial: 30 s
%       2. parallel: 10 s
%

b_pole=(abs(El)>asin(z_max_nan));      % bool to bad region around poles
nShots=size(K,1);

%%% 1. serial
% N_halo=cell(1,2);
% for mm=1:2
%     N_halo{mm}=cellfun(@(k) haloZoneCount(k,Az,El,...
%         sig_mode,lim_mode,count_mode),K(:,mm),'UniformOutput',false);
%     
%     N_halo{mm}=cat(3,N_halo{mm}{:});   % concatenate cell-array (shots) to matrix
%     N_halo{mm}(repmat(b_pole,[1,1,nShots]))=NaN;   % handle empty polar regions
% end

%%% 2. parallel
N_halo=cell(1,2);
for mm=1:2
    N_halo_temp=NaN([size(Az),nShots]);
    parfor ii=1:nShots
        N_halo_temp(:,:,ii)=haloZoneCount(K{ii,mm},Az,El,sig_mode,lim_mode,count_mode);
    end
    N_halo{mm}=N_halo_temp;
    
    % handle empty polar regions
    N_halo{mm}(repmat(b_pole,[1,1,nShots]))=NaN;   
end


%% Bell correlations
E=haloBellCorr(N_halo{:});

% report summary
nz_tot=numel(E);    % total number of momenta zones scanned
nz_nan=sum(isnan(E(:)));  % number of zones with undefined correlation
nz_valid=nz_tot-nz_nan;

E_mean=mean(E(:),'omitnan');
E_err=std(E(:),'omitnan')/sqrt(nz_valid);

fprintf('NaN fraction = %0.2g\n',nz_nan/nz_tot);
fprintf('avg E corr = %0.2g (%0.1g)\n',E_mean,E_err);


% visualisation
% E corr
figure;
plotFlatMapWrappedRad(Az,El,E,'eckert4');
cbar=colorbar('southoutside');
cbar.Label.String='E(P)';

%%% summary
% histogram of corrs
figure;
histogram(E(:),'Normalization','pdf');
xlabel('$E$');
ylabel('PDF');

%% end of code
toc(tstart);    % report elapsed time