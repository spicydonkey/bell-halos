% Load Raman data

% User config
% path_config='config_loop2_20171031_1.m';
% path_config='config_loop2_20171103_1.m';
% path_config='config_loop2_4ms_20171106_2.m';
% path_config='config_loop2_test.m';
path_config='config_loop_v2_20180119.m';

% vars to save to output
vars_save={'configs','path_config',...
    'fullpath_config',...
    'fname_save','path_save',...
    'nloop',...
    'bec_cent','Nsc_m','dk_m',...
    'nazim','nelev',...
    'Az','El',...
    'mf','mbool','nloop_m',...
    'amp_m','ver_m',...
    'nn_halo_m','P_rabi_m','th_rabi_m',...
    'ct_dth','nster_dth','dTheta',...
    };

%% Set-up script
t_main_start=tic;

% configure
dirsource=fileparts(mfilename('fullpath'));     % device independent source directory
fullpath_config=fullfile(dirsource,'config',path_config);   % full path to configs dir
run(fullpath_config);

% update configuration
configs=updateConfig(configs);

% set up this run's ID and misc paths
run_id=getdatetimestr;
fname_save=[mfilename,'__',run_id];
path_save=fullfile(fileparts(configs.files.path_loop),'arch',fname_save);
if configs.flags.savedata || configs.flags.savefigs
    mkdir(path_save);   % create dir to output any files to save to disk
end

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

%% load halos and get counts in zone
% common spherical momentum zones
nazim=configs.zone.nazim;
nelev=configs.zone.nelev;
[Az,El]=sphgrid(nazim,nelev);

% preallocate vars
mf=NaN(nloop,1);
amp=NaN(nloop,1);
ver=NaN(nloop,1);
nan_alloc=NaN(size(Az,1),size(Az,2),nloop);
nn_halo=cell(1,2);
for ii=1:2
    nn_halo{ii}=nan_alloc;
end

% TODO - need to check if BEC centers vary significantly between exps
bec_cent=cell(nloop,1);

% halo summary
Nsc=NaN(nloop,2);   % esimate of no. scattered atoms in halo - not very useful in this analysis
dk=NaN(nloop,2);    % halo rms width in k-space

% get count distribution in halo
for ii=1:nloop
    tdir=dloop{ii};
    [tmf,tamp,tver]=getLoopInfo(tdir);
    
    %% pre-process zxy-counts for this run
    dpath=fullfile(tdir,'d');
    [txy,fout]=load_txy(dpath,[],configs.load.window,...
        configs.load.mincount,configs.load.maxcount,...
        configs.load.rot_angle,1,...
        configs.flags.verbose,0);
    zxy=txy2zxy(txy);
    nshots=size(zxy,1);
    
    %% get halos for mf=0,1
    halo_k=cell(nshots,2);
    tNsc=NaN(1,2);
    tdk=NaN(1,2);
    tbec_cent=cell(1,2);
    hfig_halos=cell(1,2);
    for jj=1:2
        [halo_k(:,jj),tNsc(jj),~,tdk(jj),~,tbec_cent{jj},hfig_halos{jj}]=halo_2bec(zxy,...
            configs.bec{jj}.pos{1},configs.bec{jj}.pos{2},...
            configs.bec{jj}.Rmax,configs.bec{jj}.r_th,...
            configs.halo{jj}.dR,configs.halo{jj}.elev_max,...
            configs.flags.verbose);
    end

    %% get count distribution around the sphere 
    % NOTE
    % - combines all shots - ignore shot-to-shot noise for theta characterisation
    tnn_halo=cell(1,2);
    for mm=1:2
        tnn_halo{mm}=haloZoneCount(vertcat(halo_k{:,mm}),...
            Az,El,...
            configs.zone.sig,configs.zone.lim,...
            configs.zone.histtype);     % counts at Az, El ndgrid zones
    end
    
    % pad bad zones with NaN
    % - [x] package into a function
    % - [x] poles
    % - [ ] bright/dark spots: spontaneous halo
    %%% poles - have been removed since it's close to BEC
    el_pol_max=asin(0.7);
    el_pol=min([el_pol_max,cellfun(@(s) s.elev_max,configs.halo)]);      % min common pol limit
    tnn_halo=cellfun(@(c) pad_sphgrid_poles(El,c,el_pol,NaN),tnn_halo,'UniformOutput',false);
    
    %% store results
    mf(ii)=tmf;
    amp(ii)=tamp;
    ver(ii)=tver;
    for jj=1:2
        nn_halo{jj}(:,:,ii)=tnn_halo{jj};
    end
    
    Nsc(ii,:)=tNsc;
    dk(ii,:)=tdk;
    bec_cent{ii}=tbec_cent;
    
    % save figs output
    if configs.flags.verbose>0
        for mm=1:2
            for ll=1:numel(hfig_halos{mm})
                % get this fig
                t_hfig=hfig_halos{mm}(ll);
                figname=sprintf('fig_preproc_%0.2g_%0.2g_%d_%d_%d',tver,tamp,tmf-1,mm-1,ll);
                if configs.flags.savefigs
                    % save fig
                    saveas(t_hfig,[fullfile(path_save,figname),'.png']);
                end
            end
        end
    end
    close(hfig_halos{:});
end
clearvars txy zxy halo_k;   % clean workspace


%% tidy loop results
% sort data indices by mF
mbool{1}=(mf==0);
mbool{2}=(mf==1);
nloop_m=cellfun(@sum,mbool);

% categorise by mf
amp_m=cell(1,2);
ver_m=cell(1,2);
nn_halo_m=cell(1,2);
Nsc_m=cell(1,2);
dk_m=cell(1,2);
for ii=1:2
    amp_m{ii}=amp(mbool{ii});
    ver_m{ii}=ver(mbool{ii});
    for jj=1:2
        nn_halo_m{ii}{jj}=nn_halo{jj}(:,:,mbool{ii});
    end
    Nsc_m{ii}=Nsc(mbool{ii},:);
    dk_m{ii}=dk(mbool{ii},:);
end

% sort by Raman amp
for ii=1:2
    % get sorted index list
    [amp_m{ii},Is]=sort(amp_m{ii});
    
    % sort rest
    ver_m{ii}=ver_m{ii}(Is);
    for jj=1:2
        nn_halo_m{ii}{jj}=nn_halo_m{ii}{jj}(:,:,Is);
    end
    Nsc_m{ii}=Nsc_m{ii}(Is,:);
    dk_m{ii}=dk_m{ii}(Is,:);
end


%% evaluate Rabi population
P_rabi_m=cell(1,2);
for ii=1:2
    tnn_tot=nn_halo_m{ii}{1}+nn_halo_m{ii}{2};  % total pop in zone
    tnn_orig=nn_halo_m{ii}{ii};     % pop of original state
    P_rabi_m{ii}=tnn_orig./tnn_tot;
end


%% Rotation angle
% TODO - Rabi oscillation fitted theta - amplitude
%
% Model selection
% [ ] simple physical model
% [ ] trend fit

%%% 1. Simple formula
% NOTE: wraps rotation angle to [0,pi]
th_rabi_m=cell(size(P_rabi_m));
for ii=1:2
    if nloop_m(ii)==0
        continue
    end
    for jj=1:nloop_m(ii)
        this_P_rabi=P_rabi_m{ii}(:,:,jj);
        this_th_rabi=2*acos(sqrt(this_P_rabi));
        th_rabi_m{ii}(:,:,jj)=this_th_rabi;
    end
end

% %%% 2. Model
% %%%% fit Rabi oscillation
% aRabiFreq=NaN(size(Az,1),size(Az,2),2);
% for kk=1:2
%     this_amp=amp_m{kk};
%     if isempty(this_amp)
%         continue
%     end
%     for ii=1:size(Az,1)
%         for jj=1:size(Az,2)
%             this_pp=squeeze(P_rabi_m{kk}(ii,jj,:));
% %             if kk==1
% %                 this_pp=1-this_pp;      % for mf=0 the pop reverses
% %             end
%             aRabiFreq(ii,jj,kk)=fitRabiOsc(this_amp,this_pp);
%         end
%     end
% end
% 
% %%%% evaluate fitted theta map
% % fitted rotation angle at the data amps
% th_rabi_m=cell(size(P_rabi_m));
% for ii=1:2
%     th_rabi_m{ii}=repmat(aRabiFreq(:,:,ii),[1,1,nloop_m(ii)]);
%     for jj=1:nloop_m(ii)
%         th_rabi_m{ii}(:,:,jj)=th_rabi_m{ii}(:,:,jj)*amp_m{ii}(jj);
%     end
% end

%% 1D zonal histogram
% create histogram bin for 1D relative angles
nbin_dth=configs.zone.nbin_dth_1d;
ed_dth=linspace(-pi,pi,nbin_dth+1);
ct_dth=ed_dth(1:end-1)+0.5*diff(ed_dth);

% create gaussian filter
g1d_hsize=configs.zone.g1d_hsize;
g1d_sigma=configs.zone.g1d_sigma;
g1d_filt=gaussFilter(g1d_hsize,g1d_sigma);

nster_dth=cell(1,2);

% evaluate relative angular map for each run
dTheta=cell(1,2);
for ii=1:2
    this_nloop=nloop_m(ii);    
    
    % calculate relative rotation angle map
    % NOTE - we assume mf=0 rotates identically to mf=1. mf=0 data has too much background
    % this is a trial dth - assuming mf=0 rotates like mf=1
    dTheta{ii}=th_rabi_m{ii}-flip_bb(th_rabi_m{ii},1);
    
    nster_dth{ii}=zeros(this_nloop,nbin_dth);     % preallocate data
    for jj=1:this_nloop        
        this_dth=dTheta{ii}(:,:,jj);

        this_nster=histsphzone(this_dth,Az,El,ed_dth);
        this_nster_filt=conv(this_nster,g1d_filt,'same');
%         this_nster_filt=this_nster;

        nster_dth{ii}(jj,:)=this_nster_filt;
    end
end


%% plot
if configs.flags.graphics
    %% Flat halo density distribution    
    t_hfig=figure();
    for tmf=1:2
        tnplot=nloop_m(tmf);
        for ii=1:tnplot
            % compare count dist between mf for this run
            clf(t_hfig);
            for jj=1:2
                subplot(1,2,jj);
                plotFlatMapWrappedRad(Az,El,nn_halo_m{tmf}{jj}(:,:,ii),'eckert4');
                
                % annotate fig
                t_cb=colorbar('southoutside');
                t_cb.Label.String='density [a.u.]';

                strTitle=sprintf('[src $m_F=%d$] $m_F=%d$, $K_R=%0.2g$',tmf-1,jj-1,amp_m{tmf}(ii));
                title(strTitle);
                
            end
            drawnow
            
            % save fig
            figname=sprintf('fig_ndist_%d_%0.2g',tmf-1,amp_m{tmf}(ii));
            if configs.flags.savefigs
                saveas(t_hfig,[fullfile(path_save,figname),'.png']);
            end
        end
    end
    
    %% Population ratio - spherical distribution
    
    t_hfig=figure();
    for tmf=1:2
        tnplot=nloop_m(tmf);
        for ii=1:tnplot
            clf(t_hfig);
            
            plot_sph_surf(Az,El,P_rabi_m{tmf}(:,:,ii));
            
            % annotate fig
            t_cb=colorbar('southoutside');
            t_cb.Label.String='P_0';
            
            strTitle=sprintf('[src $m_F=%d$] $K_R=%0.2g$',tmf-1,amp_m{tmf}(ii));
            title(strTitle);
            
            drawnow
            
            % save fig
            figname=sprintf('fig_pop_%d_%0.2g',tmf-1,amp_m{tmf}(ii));
            if configs.flags.savefigs
                saveas(t_hfig,[fullfile(path_save,figname),'.png']);
            end
        end
    end
    
    %% Theta - spherical distribution
    t_hfig=figure();
    for tmf=1:2
        tnplot=nloop_m(tmf);
        
        for ii=1:tnplot
            clf(t_hfig);
            
            plot_sph_surf(Az,El,th_rabi_m{tmf}(:,:,ii));
              
            % annotate fig
            t_cb=colorbar('southoutside');
            t_cb.Label.String='\theta';
            
            strTitle=sprintf('[src $m_F=%d$] $K_R=%0.2g$',tmf-1,amp_m{tmf}(ii));
            title(strTitle);
            
            drawnow
            
            % save fig
            figname=sprintf('fig_theta_%d_%0.2g',tmf-1,amp_m{tmf}(ii));
            if configs.flags.savefigs
                saveas(t_hfig,[fullfile(path_save,figname),'.png']);
            end
        end
    end
    
    %% Rabi oscillation at selected momentum modes
    % define modes to plot Rabi process
    ndiv_az=configs.zone.ndiv_az;
    ndiv_el=configs.zone.ndiv_el;
    az_idx=1:ceil(nazim/ndiv_az):nazim-1;
    el_idx=1:ceil(nelev/ndiv_el):nelev-1;
    [Az_idx,El_idx]=ndgrid(az_idx,el_idx);
    
    % config - fig annotations
    pcolors=distinguishable_colors(ndiv_az*ndiv_el);
    pmarkers={'o','^'};
    plinestyles={'-','--'};
    
    % Population
    for tmf=1:2
        if isempty(P_rabi_m{tmf})
            continue
        end
        
        t_hfig=figure();
        for ii=1:numel(Az_idx)
            hold on;
            plot(amp_m{tmf},squeeze(P_rabi_m{tmf}(Az_idx(ii),El_idx(ii),:)),...
                'LineStyle',plinestyles{tmf},...
                'Marker',pmarkers{tmf},...
                'Color',pcolors(ii,:));
        end
        
        % annotate fig
        box on;
        hold off;
        xlabel('Raman Amplitude');
        ylabel('$P_{0}$');
        strTitle=sprintf('[src $m_F=%d$]',tmf-1);
        title(strTitle);
        
        drawnow
        
        % save fig
        figname=sprintf('fig_rabipop_%d',tmf-1);
        if configs.flags.savefigs
            saveas(t_hfig,[fullfile(path_save,figname),'.png']);
        end
    end
    
    % theta
    for tmf=1:2
        if isempty(th_rabi_m{tmf})
            continue
        end
        
        t_hfig=figure();
        for ii=1:numel(Az_idx)
            hold on;
            plot(amp_m{tmf},squeeze(th_rabi_m{tmf}(Az_idx(ii),El_idx(ii),:)),...
                'LineStyle',plinestyles{tmf},...
                'Marker',pmarkers{tmf},...
                'Color',pcolors(ii,:));
        end
        
        % annotate fig
        box on;
        hold off;
        xlabel('Raman Amplitude');
        ylabel('$\Theta$');
        strTitle=sprintf('[src $m_F=%d$]',tmf-1);
        title(strTitle);
        
        drawnow
        
        % save fig
        figname=sprintf('fig_theta_%d',tmf-1);
        if configs.flags.savefigs
            saveas(t_hfig,[fullfile(path_save,figname),'.png']);
        end
    end
    
    %% Relative rotation angle map
    % TODO
    % [ ] test
        
    t_hfig=figure();
    for tmf=1:2
        tnplot=nloop_m(tmf);
        for ii=1:tnplot
            clf(t_hfig);

            plotFlatMapWrappedRad(Az,El,dTheta{tmf}(:,:,ii),'eckert4');

            % annotate fig
            t_cb=colorbar('southoutside');
            t_cb.Label.String='\Delta\Theta';

            strTitle=sprintf('[src $m_F=%d$] $K_R=%0.2g$',tmf-1,amp_m{tmf}(ii));
            title(strTitle);
                
            drawnow
            
            % save fig
            figname=sprintf('fig_dthmap_%d_%0.2g',tmf-1,amp_m{tmf}(ii));
            if configs.flags.savefigs
                saveas(t_hfig,[fullfile(path_save,figname),'.png']);
            end
        end
    end

    %% Relative angle histogram
    for tmf=1:2
        if nloop_m(tmf)==0
            continue
        end
        
        t_hfig=figure();
        
        cc=distinguishable_colors(nloop_m(tmf));

        pp=[];
        for ii=1:nloop_m(tmf)
            ss=nster_dth{tmf}(ii,:);
            tstr=sprintf('%0.2g',amp_m{tmf}(ii));
            
            hold on;
            pp(ii)=plot(ct_dth,ss,...
                'Color',cc(ii,:),'LineWidth',1.5,...
                'DisplayName',tstr);
        end
        
        % mark CHSH angles
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
        ax=gca;
        xlabel('$\Delta\psi$');
        ylabel('Solid angle in halo [sr]');
        box on;
        xlim([-pi,pi]);
        leg=legend(pp);
        leg.Title.String='Raman amp.';
        ax.FontSize=12;
        leg.FontSize=10;
        
        drawnow
        
        % save fig
        figname=sprintf('fig_histdTheta_%d',tmf-1);
        if configs.flags.savefigs
            saveas(t_hfig,[fullfile(path_save,figname),'.png']);
        end
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
    save([path_save,'.mat'],varsExist{:});
end

%% END
t_main_end=toc(t_main_start);
disp('-----------------------------------------------');
fprintf('Total elapsed time (s): %7.1f\n',t_main_end);
disp('===================ALL TASKS COMPLETED===================');
