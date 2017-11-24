% prototyping: load Bell data

%% configure
%%% FLAGS
configs.flags.verbose=2;
configs.flags.savedata=0;       % TODO - req'd currently since each stage passes data by save/load to disk
configs.flags.archive_txy=0;        % archives loaded TXY as .mat file for future reuse
configs.flags.graphics=1;       % toggle to control graphics/plotting options
configs.flags.build_txy=1;

%%% Load
% dpath='\\amplpc29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\bell_v2\bell\run_3ms_0.38_lown_1\d';
% configs.load.id=[];         % file id numbers to use for analysis
% 
dpath='\\amplpc29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\bell_v2\bell\run_3ms_0.38_lown_2\d';
configs.load.id=[];         % file id numbers to use for analysis

configs.load.mincount=500;         % min counts in window - 0 for no min
configs.load.maxcount=1300;          % max counts in window - Inf for no max

% Detector/trap alignment
configs.load.rot_angle=0.61;

% TXY window - region of interest ( [] --> no crop )
configs.load.window{1}=[0.38,0.44];      % T [s]
configs.load.window{2}=[-35e-3,35e-3];    % X [m]
configs.load.window{3}=[-35e-3,35e-3];    % Y [m]

%%% HALO
% TODO - update with v2 loop info
%%% HALO PARAMS: BEC counts + oscillation removal for broad capture of halos
configs.bec{1}.pos{1}=[1.6525,-3.3e-3,4.5e-3];   % approx condensate locations (z,x,y)
configs.bec{1}.pos{2}=[1.7056,-2.8e-3,4.7e-3];
configs.bec{1}.Rmax=8e-3;      % max condensate sph radius
configs.bec{1}.r_th=8e-3;
configs.halo{1}.dR=[0.2,0.1];      % broad radial mask fractional width (in/out)
configs.halo{1}.zcap=0.8;   % z-cutoff (fractional wrt radius)
configs.halo{1}.string='$m_F=0$';

configs.bec{2}.pos{1}=[1.5781,-3.4e-3,3.4e-3];   % approx condensate locations (z,x,y)
configs.bec{2}.pos{2}=[1.6289,-2.9e-3,3.9e-3];
configs.bec{2}.Rmax=8e-3;      % max condensate sph radius
configs.bec{2}.r_th=8e-3;     % BEC tail radial frac diff
configs.halo{2}.dR=[0.4,0.2];      % broad radial mask fractional width (in/out)
configs.halo{2}.zcap=0.8;   % z-cutoff (fractional wrt radius)
configs.halo{2}.string='$m_F=1$';

%%% MISCELLANEOUS
configs.misc.vel_z=9.8*0.416;    % atom free-fall vert v at detector hit for T-to-Z conversion;
vz=configs.misc.vel_z;


%% main
configs=updateConfig(configs,1);

%% load txy
[txy,fout]=load_txy(dpath,configs.load.id,configs.load.window,...
    configs.load.mincount,configs.load.maxcount,...
    configs.load.rot_angle,1,...
    configs.flags.verbose,1);
zxy=txy2zxy(txy);

nShots=numel(zxy);

%% get halos for mf=0,1
halo_k=cell(nShots,2);
Nsc_mean=NaN(1,2);
Nsc_sdev=NaN(1,2);
dk_rms_mean=NaN(1,2);
dk_rms_sdev=NaN(1,2);
bec_cent=cell(1,2);
hfig_halos=cell(1,2);
% for ii=1:2
for ii=1:2
    [halo_k(:,ii),Nsc_mean(ii),Nsc_sdev(ii),dk_rms_mean(ii),dk_rms_sdev(ii),bec_cent{ii},hfig_halos{ii}]=halo_2bec(zxy,...
        configs.bec{ii}.pos{1},configs.bec{ii}.pos{2},...
        configs.bec{ii}.Rmax,configs.bec{ii}.r_th,...
        configs.halo{ii}.dR,configs.halo{ii}.elev_max,...
        configs.flags.verbose);
end