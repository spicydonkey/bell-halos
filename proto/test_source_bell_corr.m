% Testing Bell correlation distribution around source halo

% load data
load('run_source1_20171212.mat');
% loads entangled source halos in K


%% configure
% scattered zones
nAz=100;
nEl=50;

az=linspace(-pi,pi,nAz+1);
az=az(1:end-1);
el=linspace(-pi/2,pi/2,nEl);

[Az,El]=ndgrid(az,el);

% counting
sig_mode=[0.03,Inf];
lim_mode=[2,Inf];
count_mode='gauss';


%% get counts in momenta modes
nShots=size(K,1);
b_bad=(abs(El)>asin(0.7));      % polar caps culled - bad spherical region

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
figure;
plotFlatMapWrappedRad(Az,El,-E,'eckert4');
cbar=colorbar('SouthOutside');
cbar.Label.String='-E(0,0)';