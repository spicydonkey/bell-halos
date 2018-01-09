% Testing Bell correlation distribution in SOURCE
%

tstart=tic;

%% load data
load('src1_data_digest.mat');
% loads entangled source halos (halo_k0)
K=halo_k0;


%% configure
% scattered modes
nAz=200;        % num of azimuthal divisions in [-pi,pi)
nEl=100;         % num of elevation divisions in [-pi/2,pi/2]
% NOTE: nAz should be EVEN to ensure flip_bb.m correctly matches
% spherically opposite modes

% counting
count_mode='gauss';
sig_mode=[0.018,Inf];
lim_mode=[3,Inf];


%% set-up
% scattered modes
az=linspace(-pi,pi,nAz+1);
az=az(1:end-1);
el=linspace(-pi/2,pi/2,nEl);

[Az,El]=ndgrid(az,el);


%% get counts in momenta modes
z_max_nan=0.7;      % max-z used to filter halos as bad
b_bad=(abs(El)>asin(z_max_nan));      % bool to bad region around poles

nShots=size(K,1);

N_halo=cell(1,2);
for mm=1:2
    N_halo{mm}=cellfun(@(k) haloZoneCount(k,nAz+1,nEl,...
        sig_mode,lim_mode,count_mode),K(:,mm),'UniformOutput',false);
    % gaussian weighted counting
    
    N_halo{mm}=cat(3,N_halo{mm}{:});   % concatenate cell-array (shots) to matrix
    
    N_halo{mm}(repmat(b_bad,[1,1,nShots]))=NaN;   % handle empty polar regions
end


%% Bell correlations
E=haloBellCorr(N_halo{:});


% report summary
nz_tot=numel(E);    % total number of momenta zones scanned
nz_nan=sum(sum(isnan(E)));  % number of zones with undefined correlation

E_mean=sum(sum(E,'omitnan'),'omitnan')/(nz_tot-nz_nan);

fprintf('NaN fraction = %0.2g\n',nz_nan/nz_tot);
fprintf('avg E corr = %0.2g\n',E_mean);


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