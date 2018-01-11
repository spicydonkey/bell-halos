% Testing Bell correlation distribution in BELL
%

tstart=tic;

%% load data
%%% Bell test
load('run2_20171213_0.mat');
% loads g2-optimised bell-halos (K_g2_opt)
K=K_g2_opt;

% % pre-radial transform
% load('run2_Knew_g2_opt_20171212.mat');
% % K is original

%%% Calibrated relative rotation - 3ms,0.382
S_loop=load('dpsi_20171213_0.mat');


%% configure
% scattered modes
nAz=100;        % num of azimuthal divisions in [-pi,pi)
nEl=50;         % num of elevation divisions in [-pi/2,pi/2]
% NOTE: nAz should be EVEN to ensure flip_bb.m correctly matches
% spherically opposite modes

% atom counting
% count_mode='gauss';
% sig_mode=[0.0133,Inf];
% lim_mode=[4.5,Inf];

count_mode='simple';
sig_mode=0.05;
lim_mode=[];

z_max_nan=0.7;      % max-z used to filter halos as bad


%% set-up
[Az,El]=sphgrid(nAz,nEl);


%% get counts in momenta modes

b_pole=(abs(El)>asin(z_max_nan));      % bool to bad region around poles

nShots=size(K,1);

N_halo=cell(1,2);
for mm=1:2
    N_halo{mm}=cellfun(@(k) haloZoneCount(k,Az,El,...
        sig_mode,lim_mode,count_mode),K(:,mm),'UniformOutput',false);
    
    N_halo{mm}=cat(3,N_halo{mm}{:});   % concatenate cell-array (shots) to matrix
    N_halo{mm}(repmat(b_pole,[1,1,nShots]))=NaN;   % handle empty polar regions
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

%%%% visualisation
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


%% Bell correlation vs. rel rotation
% generate interpolated rotation map
dPsi_interp=interpn(S_loop.Az,S_loop.El,S_loop.dPsi,Az,El);

% visualise
figure;
plotFlatMapWrappedRad(Az,El,dPsi_interp,'eckert4');
cbar=colorbar('southoutside');
cbar.Label.String='\Delta\psi [rad]';

% clean data
EE=E(:);
DD=dPsi_interp(:);
DD=wrapToPi(DD);
[DD,Is]=sort(DD,'ascend');      % sort by rel rot angle
EE=EE(Is);

% plot: E vs dpsi
figure;
plot(DD,EE,'b.');
xlim([-pi,pi]);
ylim([-1,1]);
xlabel('$\Delta\psi$');
ylabel('$E$');


%% smoothing
n_dpsi_bin=31;
dpsi_ed=linspace(-pi,pi,n_dpsi_bin+1);
dpsi_ct=dpsi_ed(1:end-1)+0.5*diff(dpsi_ed);

DD_sm=zeros(2,n_dpsi_bin);   % mean and serr
EE_sm=zeros(2,n_dpsi_bin);  

for ii=1:n_dpsi_bin
    bool_bin=(DD>dpsi_ed(ii)&DD<dpsi_ed(ii+1));     % boolean for items in bin
    n_items=sum(bool_bin);      % number of items in this bin
    
    % get values
    d_bin=DD(bool_bin);
    e_bin=EE(bool_bin);
    
    % statistics
    DD_sm(1,ii)=mean(d_bin);
    DD_sm(2,ii)=std(d_bin)/sqrt(n_items);
    
    EE_sm(1,ii)=mean(e_bin,'omitnan');
    EE_sm(2,ii)=std(e_bin,'omitnan')/sqrt(n_items);
end

%%% plot
figure;
% plot(DD_sm(1,:),EE_sm(1,:),'b.');
hh=ploterr(DD_sm(1,:),EE_sm(1,:),DD_sm(2,:),EE_sm(2,:),'o','hhxy',0);
set(hh(1),'Color','b','MarkerFaceColor','w');	% DATAPOINT
set(hh(2),'Color','b','DisplayName','');                  % Y-err
set(hh(3),'Color','b','DisplayName','');                  % X-err

% THEORIES
% LHV
hold on;
hLHV=line([-pi,0,pi],[1,-1,1],'LineStyle','--','Color',[0.5,0.5,0.5],'DisplayName','LHV');
hold off;

% QM
hold on;
dth_qm=linspace(-pi,pi,100);
E_qm=-cos(dth_qm);
hQM=plot(dth_qm,E_qm,'LineStyle','-','Color','k','DisplayName','QM');
hold off;

xlim([-pi,pi]);
ylim([-1,1]);
xlabel('$\Delta\psi$');
ylabel('$E$');

%% end of code
toc(tstart);    % report elapsed time