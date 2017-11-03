% Load Raman data

% User config
% path_config='config_loop2_20171031_1.m';
path_config='config_loop2_20171103_1.m';


% vars to save to output
vars_save={'configs','path_config',...
    'nloop',...
    'nazim','nelev',...
    'Az','El',...
    'mf','mbool','nloop_m',...
    'amp_m','ver_m',...
    'P_rabi_m','th_rabi_m',...
    'ct_dth','ster_dth','DTHETA',...
    };

%% Main
t_main_start=tic;

% configure
dirsource=fileparts(mfilename('fullpath'));     % device independent source directory
fullpath_config=fullfile(dirsource,'config',path_config);   % full path to configs dir
run(fullpath_config);

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
nan_alloc=NaN(size(Az,1),size(Az,2),nloop);
nn_halo=cell(1,2);
for ii=1:2
    nn_halo{ii}=nan_alloc;
end
P_rabi=nan_alloc;
th_rabi=nan_alloc;

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
    
    % analyse the Rabi state flopping for this parameter set (pop is for
    % mf=1)
    [~,~,tP_rabi,tnn_halo]=rabiAnalyse(halo_k,configs.zone.nazim,configs.zone.nelev,...
        configs.zone.sig,configs.zone.lim,...
        configs.zone.histtype);
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
    for jj=1:2
        nn_halo{jj}(:,:,ii)=tnn_halo{jj};
    end
    P_rabi(:,:,ii)=tP_rabi;
    th_rabi(:,:,ii)=tth_rabi;
    
    bec_cent{ii}=tbec_cent;
end
mbool{1}=(mf==0);
mbool{2}=(mf==1);
nloop_m=cellfun(@sum,mbool);

%% clean loop results
% categorise to mf
amp_m=cell(1,2);
ver_m=cell(1,2);
nn_halo_m=cell(1,2);
P_rabi_m=cell(1,2);
th_rabi_m=cell(1,2);

for ii=1:2
    amp_m{ii}=amp(mbool{ii});
    ver_m{ii}=ver(mbool{ii});
    for jj=1:2
        nn_halo_m{ii}{jj}=nn_halo{jj}(:,:,mbool{ii});
    end
    P_rabi_m{ii}=P_rabi(:,:,mbool{ii});
    th_rabi_m{ii}=th_rabi(:,:,mbool{ii});
end

% sort by Raman amp
for ii=1:2
    [amp_m{ii},Is]=sort(amp_m{ii});
    ver_m{ii}=ver_m{ii}(Is);
    for jj=1:2
        nn_halo_m{ii}{jj}=nn_halo_m{ii}{jj}(:,:,Is);
    end
    P_rabi_m{ii}=P_rabi_m{ii}(:,:,Is);
    th_rabi_m{ii}=th_rabi_m{ii}(:,:,Is);
end

%% 1D zonal histogram
nzones=numel(Az);
sterPerZone=(4*pi)/nzones;

% create histogram bin for 1D relative angles
nbin_dth=configs.zone.nbin_dth_1d;
ed_dth=linspace(-pi,pi,nbin_dth);
ct_dth=ed_dth(1:end-1)+0.5*diff(ed_dth);

% create gaussian filter
g1d_hsize=configs.zone.g1d_hsize;
g1d_sigma=configs.zone.g1d_sigma;
g1d_filt=gaussFilter(g1d_hsize,g1d_sigma);

ster_dth=cell(1,2);
DTHETA=cell(1,2);

for ii=1:2
    this_nloop=nloop_m(ii);
    % preallocate 
    ster_dth{ii}=zeros(this_nloop,length(ct_dth));
    
    % TODO 
    % this is a trial dth - assuming mf=0 rotates like mf=1
    DTHETA{ii}=th_rabi_m{ii}-flip_bb(th_rabi_m{ii});
%     DTHETA{ii}=NaN(size(th_rabi_m{ii}));
%     for jj=1:this_nloop
%         DTHETA{ii}(:,:,jj)=th_rabi_m{ii}(:,:,jj)-flip_bb(th_rabi_m{ii}(:,:,jj));
%     end
    
    for jj=1:this_nloop        
        this_dth=DTHETA{ii}(:,:,jj);
        this_dth=this_dth(:);
        
        tnn=histcounts(this_dth,ed_dth);
        tster=tnn*sterPerZone;
        tster=conv(tster,g1d_filt,'same');
        ster_dth{ii}(jj,:)=tster;
    end
end

%% plot
if configs.flags.graphics
    %% Flat halo density distribution
    % TODO
    % annotate plots 
    
    for tmf=1:2
        tnplot=nloop_m(tmf);
        
        for ii=1:tnplot
            % compare count dist between mf for this run
            figure();
            for jj=1:2
                subplot(2,1,jj);
                plotFlatMap(rad2deg(El),rad2deg(Az),nn_halo_m{tmf}{jj}(:,:,ii),'eckert4');
            end
        end
    end
    
    %% Spherical distribution
    % Population ratio
    figure();
    tmf=1;
    tnplot=nloop_m(tmf);
    for ii=1:tnplot
        subplot(1,tnplot,ii);
        plot_sph_surf(Az,El,P_rabi_m{tmf}(:,:,ii));
        colorbar('southoutside');
    end
    
    figure();
    tmf=2;
    tnplot=nloop_m(tmf);
    for ii=1:tnplot
        subplot(1,tnplot,ii);
        plot_sph_surf(Az,El,P_rabi_m{tmf}(:,:,ii));
        colorbar('southoutside');
    end
    
    % Theta
    figure();
    tmf=1;
    tnplot=nloop_m(tmf);
    for ii=1:tnplot
        subplot(1,tnplot,ii);
        plot_sph_surf(Az,El,th_rabi_m{tmf}(:,:,ii));
        colorbar('southoutside');
    end

    figure();
    tmf=2;
    tnplot=nloop_m(tmf);
    for ii=1:tnplot
        subplot(1,tnplot,ii);
        plot_sph_surf(Az,El,th_rabi_m{tmf}(:,:,ii));
        colorbar('southoutside');
    end
    
    %% Modes
    % define modes to plot Rabi process
    ndiv_az=configs.zone.ndiv_az;
    ndiv_el=configs.zone.ndiv_el;
    az_idx=1:ceil(nazim/ndiv_az):nazim-1;
    el_idx=1:ceil(nelev/ndiv_el):nelev-1;
    [Az_idx,El_idx]=ndgrid(az_idx,el_idx);
    
    % annotation config
    pcolors=distinguishable_colors(ndiv_az*ndiv_el);
    pmarkers={'o','^'};
    plinestyles={'-','--'};
    
    % Population
    figure(); clf;
    for ii=1:numel(Az_idx)
        for jj=1:2
            if isempty(P_rabi_m{jj})
                continue
            end
            hold on;
            plot(amp_m{jj},squeeze(P_rabi_m{jj}(Az_idx(ii),El_idx(ii),:)),...
                'LineStyle',plinestyles{jj},...
                'Marker',pmarkers{jj},...
                'Color',pcolors(ii,:));
        end
    end
    box on;
    hold off;
    xlabel('Raman Amplitude');
    ylabel('$P(\uparrow)$');
    
    % theta
    figure(); clf;
    for ii=1:numel(Az_idx)
        for jj=1:2
            if isempty(th_rabi_m{jj})
                continue
            end
            hold on;
            plot(amp_m{jj},squeeze(th_rabi_m{jj}(Az_idx(ii),El_idx(ii),:)),...
                'LineStyle',plinestyles{jj},...
                'Marker',pmarkers{jj},...
                'Color',pcolors(ii,:));
        end
    end
    box on;
    hold off;
    xlabel('Raman Amplitude');
    ylabel('$\Theta$');
    
    %% Relative angle histogram
    for mm=1:2
        cc=distinguishable_colors(nloop_m(mm));
        
        h_ss(mm)=figure(); hold on;
        
        pp=[];
        for ii=1:nloop_m(mm)
            ss=ster_dth{mm}(ii,:);
            
            figure(h_ss(mm));
            tstr=sprintf('%0.2g',amp_m{mm}(ii));
            pp(ii)=plot(ct_dth,ss,...
                'Color',cc(ii,:),'LineWidth',1.5,...
                'DisplayName',tstr);
        end
        figure(h_ss(mm)); hold off;
        
        % test angles - lines
        testAngles=[-pi/4,pi/4,3/4*pi];
        
        hold on;
        yylim=ylim;
        ylim(yylim);        % fix ylim
        plines=[];
        gray_col=0.8*ones(1,3);         % gray data points
        for ii=1:numel(testAngles)
            plines(ii)=line(testAngles(ii)*[1,1],yylim,...
                'Color',gray_col,'LineStyle','--','LineWidth',2);
            uistack(plines(ii),'bottom');
        end
        
        % annotation
        figure(h_ss(mm));
        ax=gca;
        xlabel('$\Delta\psi$');
        ylabel('Solid angle in halo [sr]');
        box on;
        xlim([-pi,pi]);
        leg=legend(pp);
        leg.Title.String='Raman amp.';
        ax.FontSize=12;
        leg.FontSize=10;
    end

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