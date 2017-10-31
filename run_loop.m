% Load Raman data

% User config
path_config='C:\Users\HE BEC\Documents\MATLAB\bell-halos\config\config_loop_20171031_1.m';

% vars to save to output
vars_save={'configs','path_config',...
    'nloop',...
    'nazim','nelev',...
    'mf','amp','ver','mbool',...
    'Az','El',...
    'P_rabi','th_rabi',...
    };

%% Main
t_main_start=tic;

% configure
run(path_config);

% set up this run's ID and misc paths
run_id=getdatetimestr;
path_save=fullfile(fileparts(configs.files.path_loop),'arch',[mfilename,'__',run_id,'.mat']);

% parse main data directory
cd(configs.files.path_loop);
dlist=dir;
nfiles=length(dlist);

% get loop runs
dloop=cell(nfiles,1);
nloop=0;
for ii=1:nfiles
    if ~isnan(getLoopInfo(dlist(ii).name))
        nloop=nloop+1;
        dloop{nloop}=dlist(ii).name;
    end
end
dloop=dloop(1:nloop);

%% process loop
% common spherical momentum zones
nazim=configs.zone.nazim;
nelev=configs.zone.nelev;
az=linspace(-pi,pi,nazim);
az=az(1:end-1);             % unique angles only
el=linspace(-pi/2,pi/2,nelev);
[Az,El]=ndgrid(az,el);

% preallocate vars
mf=NaN(nloop,1);
amp=NaN(nloop,1);
ver=NaN(nloop,1);
P_rabi=NaN(size(Az,1),size(Az,2),nloop);
th_rabi=NaN(size(Az,1),size(Az,2),nloop);

% TODO - need to check if BEC centers vary significantly between exps
bec_cent=cell(nloop,1);

for ii=1:nloop
    tdir=dloop{ii};
    [tmf,tamp,tver]=getLoopInfo(tdir);
    
    % pre-process zxy-counts for this run
    dpath=fullfile(tdir,'d');
    [txy,fout]=load_txy(dpath,[],configs.load.window,...
        configs.load.mincount,configs.load.maxcount,...
        configs.load.rot_angle,1,...
        configs.flags.verbose,0);
    zxy=txy2zxy(txy);
    clearvars txy;
    nshots=size(zxy,1);
    
%     figure(1);
%     plot_zxy(zxy);
    
    % get halos for mf=0,1
    halo_k=cell(nshots,2);
    tbec_cent=cell(1,2);
    for jj=1:2
        [halo_k(:,jj),tbec_cent{jj}]=halo_2bec(zxy,...
            configs.bec{jj}.pos{1},configs.bec{jj}.pos{2},...
            configs.bec{jj}.Rmax,configs.bec{jj}.dR_tail,...
            configs.halo{jj}.dR,configs.halo{jj}.zcap,...
            configs.flags.verbose);
    end
    clearvars zxy;
    
%     figure();
%     plot_zxy(halo_k);
%     axis equal;
    
    % analyse the Rabi state flopping for this parameter set (pop is for
    % mf=1)
    tP_rabi=rabiAnalyse(halo_k,configs.zone.nazim,configs.zone.nelev,...
        configs.zone.sig,configs.zone.lim);
    clearvars halo_k;
    
    % state rotation angle at each momentum mode
    % formula: Pr(ORIGIGNAL_STATE)=cos(ROT_ANGLE/2)^2
    %   ROT_ANGLE >= 0 by convention
    % TODO - can only tell rotation angle modulo pi
    tth_rabi=2*acos(sqrt(tP_rabi));
    if tmf==0
        tth_rabi=pi-tth_rabi;
    end
    
    % store results
    mf(ii)=tmf;
    amp(ii)=tamp;
    ver(ii)=tver;
    P_rabi(:,:,ii)=tP_rabi;
    th_rabi(:,:,ii)=tth_rabi;
    
    bec_cent{ii}=tbec_cent;
end
mbool{1}=(mf==0);
mbool{2}=(mf==1);


%% process loop results
% categorise to mf
amp_m=cell(1,2);
ver_m=cell(1,2);
P_rabi_m=cell(1,2);
th_rabi_m=cell(1,2);

for ii=1:2
    amp_m{ii}=amp(mbool{ii});
    ver_m{ii}=ver(mbool{ii});
    P_rabi_m{ii}=P_rabi(:,:,mbool{ii});
    th_rabi_m{ii}=th_rabi(:,:,mbool{ii});
end

% sort by Raman amp
for ii=1:2
    [amp_m{ii},Is]=sort(amp_m{ii});
    ver_m{ii}=ver_m{ii}(Is);
    P_rabi_m{ii}=P_rabi_m{ii}(:,:,Is);
    th_rabi_m{ii}=th_rabi_m{ii}(:,:,Is);
end

%% plot
if configs.flags.graphics
    % Spherical distribution
    figure(1);
    for ii=1:nloop
        subplot(1,nloop,ii);
        plot_sph_surf(Az,El,P_rabi(:,:,ii));
        colorbar('southoutside');
    end
    
    figure(2);
    for ii=1:nloop
        subplot(1,nloop,ii);
        plot_sph_surf(Az,El,th_rabi(:,:,ii));
        colorbar('southoutside');
    end
    
    % Modes
    % define modes to plot Rabi process
    ndiv_az=configs.zone.ndiv_az;
    ndiv_el=configs.zone.ndiv_el;
    az_idx=1:ceil(nazim/ndiv_az):nazim;
    el_idx=1:ceil(nelev/ndiv_el):nelev;
    [Az_idx,El_idx]=ndgrid(az_idx,el_idx);
    
    % annotation config
    pcolors=distinguishable_colors(ndiv_az*ndiv_el);
    pmarkers={'o','^'};
    plinestyles={'-','--'};
    
    % Population
    figure(3); clf;
    for ii=1:numel(Az_idx)
        for jj=1:2
            hold on;
            %         plot(amp(mbool{jj}),squeeze(P_rabi(Az_idx(ii),El_idx(ii),mbool{jj})),'o--');
            plot(amp(mbool{jj}),squeeze(P_rabi(Az_idx(ii),El_idx(ii),mbool{jj})),...
                'LineStyle',plinestyles{jj},...
                'Marker',pmarkers{jj},...
                'Color',pcolors(ii,:));
        end
    end
    box on;
    hold off;
    xlabel('Raman Amplitude');
    ylabel('$P(\uparrow)$');
    
    % mF rotation angle
    figure(4); clf;
    for ii=1:numel(Az_idx)
        for jj=1:2
            hold on;
            %         plot(amp(mbool{jj}),squeeze(th_rabi(Az_idx(ii),El_idx(ii),mbool{jj})),'o--');
            plot(amp(mbool{jj}),squeeze(th_rabi(Az_idx(ii),El_idx(ii),mbool{jj})),...
                'LineStyle',plinestyles{jj},...
                'Marker',pmarkers{jj},...
                'Color',pcolors(ii,:));
        end
    end
    box on;
    hold off;
    xlabel('Raman Amplitude');
    ylabel('$\Theta$');
    
end

%% Save results
if configs.flags.savedata
    % TODO - package this into a function
    % TODO - if datafile already exists, MOVE it
    
    % parse list of vars user wants to save and check against workspace if exists
    varsExist=cell(size(vars_save));
    varCounter=0;
    for i = 1:length(vars_save)
        tVar=vars_save{i};
        if ~exist(tVar,'var')
            warning(['Variable "',tVar,'" does not exist.']);
            continue;
        end
        varCounter=varCounter+1;
        varsExist{varCounter}=tVar;
    end
    varsExist=varsExist(1:varCounter);
    
    dir_save=fileparts(path_save);
    if ~exist(dir_save,'dir')
        warning('Directory to save data %s does not exist. Creating directory.',dir_save);
        mkdir(dir_save);
    end
    save(path_save,varsExist{:});
end

%% END
t_main_end=toc(t_main_start);
disp('-----------------------------------------------');
fprintf('Total elapsed time (s): %7.1f\n',t_main_end);
disp('===================ALL TASKS COMPLETED===================');