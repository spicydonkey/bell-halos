% Get mode occupancy of halos generated from the Bell experiment
%

% config
txy_window={[0.38,0.44],[-35,35]*1e-3,[-35,35]*1e-3};

bec_pos_1=[1.6543,-3.5e-3,4.0e-3];
bec_pos_2=[1.7067,-3.0e-3,4.9e-3];
bec_r=8e-3;
r_th=bec_r;
halo_dR=0.2;
elev_max=asin(0.8);

dir_data='C:\Users\David\hebec\bell\bell_v2\loop\1_0.382_v2.1';
nShotsLatest=50;

halo_dk_rms=0.033;      % halo width rms
wbb=0.027;         % rms mom-width of source BEC. sigk ~ wbb/1.1


% get all reconstructed TXY shot IDs from data directory
files_txy=dir(fullfile(dir_data,'d_txy_forc*.txt'));
dfiles_id=sort(cellfun(@(x) dfile2id(x),{files_txy(:).name}'));  % get sorted shotID array

% get recent shots of interest
dfiles_id=dfiles_id(1:end-1);   % ignore the latest shot - it may be being built
% NOTE - we just take the defined number of successfully reconstructed
% shots - so it's your job to make sure AutoConvert is happy
ll=length(dfiles_id);
id_load=dfiles_id(max(1,ll-nShotsLatest):end);

% load txy from latest shots
zxy=txy2zxy(load_txy(fullfile(dir_data,'d'),id_load,txy_window,...
    0,Inf,0.61,0,0,0));
% zxy=txy2zxy(txy);

% collate all shots
% NOTE - we are ignoring shot-to-shot oscillation
%zxy=vertcat(zxy{:});

% capture mF=0 halo with polar marker BECs
%   conveniently returns scattered number in halo too!
[~,Nsc_avg,Nsc_std]=halo_2bec(zxy,bec_pos_1,bec_pos_2,bec_r,r_th,halo_dR,elev_max,0);
Nsc_unc=Nsc_std/Nsc_avg;

% evaluate halo mode occupancy
mocc_halo=halo_mocc(1,halo_dk_rms,Nsc_avg,wbb/1.1);

% report estimated mode occupancy
fprintf('halo (mf=0) mode occupancy = %0.2g (%0.1g)\n',mocc_halo,Nsc_unc*mocc_halo);
