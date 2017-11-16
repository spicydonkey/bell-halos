% Get mode occupancy of halos generated from the Bell experiment
%

% config
txy_window={[0.38,0.44],[-35,35]*1e-3,[-35,35]*1e-3};

bec_pos_1=[1.5781,-3.1e-3,2.8e-3];
bec_pos_2=[1.6297,-3.2e-3,3.6e-3];
bec_r=8e-3;
r_th=bec_r;
halo_dR=0.2;
elev_max=asin(0.8);

dir_data='/path/to/txy/';
nShotsLatest=10;

halo_dk_rms=0.033;      % halo width rms
bec_sigk=0.027;         % rms mom-width of source BEC. sigk ~ wbb/1.1


% get all reconstructed TXY shot IDs from data directory
files_all=dir(dir_data);
% dfile2id

% load txy from latest shots
% load_txy

% collate all shots

% capture mF=0 halo with polar marker BECs
%   conveniently returns scattered number in halo too!
% halo_2bec

% evaluate halo mode occupancy
% halo_mocc

% report estimated mode occupancy
