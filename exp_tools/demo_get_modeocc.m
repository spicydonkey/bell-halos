% Get mode occupancy of halos generated from the Bell experiment
%

% USER CONFIG
verbose=1;

txy_window={[0.38,0.44],[-35,35]*1e-3,[-35,35]*1e-3};

bec_pos_1=[1.6524,-3.3e-3,4.6e-3];
bec_pos_2=[1.7055,-2.8e-3,4.8e-3];
bec_r=8e-3;
r_th=bec_r;
halo_dR=0.1;        % mf=0 captures very well
elev_max=asin(0.8);

dir_data='\\amplpc29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output';
nShotsLatest=50;

halo_dk_rms=0.033;      % halo width rms
wbb=0.03;         % rms mom-width of source BEC. sigk ~ wbb/1.1
% data from dist halo


% close figs
close all;

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
[~,Nsc_avg,Nsc_std]=halo_2bec(zxy,bec_pos_1,bec_pos_2,bec_r,r_th,halo_dR,elev_max,verbose);
Nsc_unc=Nsc_std/Nsc_avg;

% evaluate halo mode occupancy
mocc_halo=halo_mocc(1,halo_dk_rms,Nsc_avg,wbb/1.1);

%% summary
fprintf('-----------------------------------------------------------------------------\n');
disp(datetime);     % tag current date-time

% analysed data
fprintf('* TXY loaded (shot ID): [%d,%d]\n',id_load(1),id_load(end));
% report estimated mode occupancy
fprintf('* halo (mf=0) mode occupancy = %0.2g (%0.1g)\n',mocc_halo,Nsc_unc*mocc_halo);

%%% end of code