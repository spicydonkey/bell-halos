%Analyse experiment on mF asymmetry
%


%% CONFIG

dbg=1;

% base directory
dir_base='X:\expdata\spinmom_bell\loop_tscan';

addpath(genpath(dir_base));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEGIN PREPROCESSING


%% 1. Load raw data
% create txy and save mat file

% or LOAD IT!
load('txy_20180129.mat');
shot_id=fout.id_ok;


%% 1.1. Load parameter log
% see load_params.m
%   params
%   id_in_params

load('param_20180131.mat');

% build boolean array selector for scanned parameter-set
b_paramset=cellfun(@(idx) ismember(shot_id,idx),...
    id_in_param,'UniformOutput',false);
b_paramset=horzcat(b_paramset{:});


%% 2. Preprocess data
%   * distinguish mf
%   * centered halo in detector coord

path_config='config_preproc.m';
run(path_config);

zxy=txy2zxy(txy);
n_mf=numel(configs.mf);

%%% BEC locations
% shot averaged - no shot-to-shot oscillation
cent_bec_avg=cell(n_mf,2);
zxy_samp=vertcat(zxy{:});
for ii=1:n_mf
    for jj=1:2
        c0=configs.mf{ii}.bec(jj,:);    % estimate
        cent_bec_avg{ii,jj}=capture_bec(zxy_samp,c0,configs.r_bec_capt);
    end
end

if dbg
    figure;
    scatter_zxy(vertcat(cent_bec_avg{:}),100,'r');
    axis equal;
    xlim(30e-3*[-1,1]);
    ylim(30e-3*[-1,1]);
end

% TODO
%   try analysis without accounting for shot oscillation
%
% % locate BECs for each shot
% cent_bec=cell(n_mf,2);
% for ii=1:n_mf
%     for jj=1:2
%         c0=configs.mf{ii}.bec(jj,:);
%         cent_bec{ii,jj}=cellfun(@(x) capture_bec(x,c0,configs.r_bec_capt),...
%             zxy,'UniformOutput',false);
%         cent_bec{ii,jj}=vertcat(cent_bec{ii,jj}{:});
%     end
% end
% % plot summary
% figure;
% for ii=1:numel(cent_bec)
%     hold on;
%     scatter_zxy(cent_bec{ii},20,'k');
% end
% axis equal;
% xlim(30e-3*[-1,1]);
% ylim(30e-3*[-1,1]);


%%% distinguish mF-species
n_shot=size(zxy,1);
zxy_mf=cell(n_shot,n_mf);
for ii=1:n_mf
    % bound Z by BECs
    box_lim={[cent_bec_avg{ii,1}(1),cent_bec_avg{ii,2}(1)],[],[]};
    zxy_mf(:,ii)=cellfun(@(x) boxcull(x,box_lim),zxy,'UniformOutput',false);
end


%%% get halo center coords
cent_halo=zeros(n_mf,3);
for ii=1:n_mf
    % estimate by mid-point of marker BECs
    cent_halo(ii,:)=mean(vertcat(cent_bec_avg{ii,:}));
end


%%% re-center to halo
zxy0_mf=cell(size(zxy_mf));
cent_bec_avg0=cell(size(cent_bec_avg));
for ii=1:n_mf
    zxy0_mf(:,ii)=boost_zxy(zxy_mf(:,ii),-cent_halo(ii,:));
    
    % BEC markers in halo-centered coord
    cent_bec_avg0(ii,:)=boost_zxy(cent_bec_avg(ii,:),-cent_halo(ii,:)); 
end

if dbg
    figure;
    colors={'r','b','k'};
    for ii=1:n_mf
        hold on;
        plot_zxy(zxy0_mf(:,ii),1e4,10,colors{ii});
        scatter_zxy(vertcat(cent_bec_avg0{ii,:}),1000,colors{ii});
    end
    axis equal;
end


%%% Clean data
zxy_halo=zxy0_mf;       % initialise halo data

% A) remove BEC + thermal atoms localised around the BECs
r_thermal=configs.r_thermal;
for ii=1:n_mf
    for jj=1:2
        c_bec=cent_bec_avg0{ii,jj};    % BEC centroid
        zxy_halo(:,ii)=cellfun(@(x) cropBall(x,r_thermal,c_bec,false),...
            zxy_halo(:,ii),'UniformOutput',false);
    end
end

if dbg
    figure;
    colors={'r','b','k'};
    for ii=1:n_mf
        hold on;
        plot_zxy(zxy_halo(:,ii),1e4,10,colors{ii});
    end
    axis equal;
end

% B) radial
r_crop=configs.r_crop;      % radial window (normalised to halo radius)
for ii=1:n_mf
    % scale to real size of halo
    tr_halo=vnorm(cent_bec_avg0{ii,1});     % BECs define poles of sphere
    r_lim=r_crop*tr_halo;
    
    zxy_halo(:,ii)=cfilter_norm(zxy_halo(:,ii),r_lim(1),r_lim(2));
end

if dbg
    figure;
    colors={'r','b','k'};
    for ii=1:n_mf
        hold on;
        plot_zxy(zxy_halo(:,ii),1e4,10,colors{ii});
    end
    axis equal;
end


%% 3. k-space
%   Unit-spherise
%

%%% ellipsoid fit
%   a fast first-attempt to spherise
k_halo=cell(size(zxy_halo));
for ii=1:n_mf
    k_halo(:,ii)=map2usph(zxy_halo(:,ii));
end

if dbg
    scatter_halo(k_halo);
end 

% %%% A) real-halo local ellipticity
% % config
% nAz=30;
% nEl=15;
% dth_cone=(2*pi/nAz);
% 
% % get spherical matp of local halo radius and thickness
% [r0,dr,az,el]=summary_disthalo_ellip(zxy_halo,[nAz,nEl],dth_cone,dbg);


%%% clean
r_lim=[0.9,1.1];
k_halo=cfilter_norm(k_halo,r_lim(1),r_lim(2));

if dbg
    scatter_halo(k_halo);
end


%% 4. Classify halo data by exp-params
% classify data in param-space
k_par=cell(n_par_iter);
for ii=1:nparam
    this_idx=num2cell(par_tab(ii,:));     % location in table 
    %cell-formated to distribute as argument list
    
    k_par{this_idx{:}}=k_halo(b_paramset(:,ii),:);      % get all halo data
    %from this param-set and store
end

if dbg
%     for ii=1:numel(k_par)
%         this_parvec=parvec_tab{ii};
%         h=summary_disthalo_ndist(k_par{ii});
%         
%         
%     end
end


%%% END OF PREPROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEGIN ANALYSIS 

path_config='config_anal.m';
run(path_config);


%% 1. Atom number distribution
% create spherical grid zones
[az,el]=sphgrid(configs.anal.nZone(1),configs.anal.nZone(2));

%%% count points around different parts of sphere
nk_sph=cell(size(k_par));
parfor ii=1:nparam
    tnk_sph=cell(n_mf,1);
    for jj=1:n_mf
        tk=vertcat(k_par{ii}{:,jj});    % shot-collated halo
        tnk_sph{jj}=haloZoneCount(tk,az,el,...
            configs.anal.sigZone,configs.anal.limZone,...
            configs.anal.typeZone);     % count points around sphere
    end
    nk_sph{ii}=cat(3,tnk_sph{:});       % sph-grid is 2D for now
    % NOTE: concatenating dim needs to be generalised to resolve radially
end

if dbg
    for ii=1:nparam
        this_parvec=parvec_tab{ii};   % get param-vec: [mf, traman]
        
        %%% plot
        figname=char(strjoin(string(this_parvec),'_'));
        figname=['N_',figname];
        
        figure('Name',figname);
        for jj=1:n_mf
            subplot(n_mf,1,jj);
            plotFlatMapWrappedRad(az,el,nk_sph{ii}(:,:,jj));
            cbar=colorbar('eastoutside');
            cbar.Label.String='N';
        end
        
        % save
        if exist('SAVEFIG','var') && SAVEFIG
            fname=fullfile(dir_base,'out','img',figname);
            print(fname,'-dpng');
        end
        
    end
end


%% 2. Population ratio and Rabi cycle

% evaluate excited state pop ratio
P=cell(size(nk_sph));       % population ratio of excited state
parfor ii=1:nparam
    this_parvec=parvec_tab{ii};   % get param-vec: [mf, traman]
    t_mf=this_parvec(1);
    
    %   mf  -->     index
    %   1   -->     1
    %   0   -->     2
    %   -1  -->     3
    idx_mf=2-t_mf;
    
    nk=nk_sph{ii};
    nk=nk(:,:,1:2);     % NOTE: ignore mf=-1 atoms!
    
    N_all=sum(nk,3);
    N_exc=N_all-nk(:,:,idx_mf);
    
    P{ii}=N_exc./N_all;
    
    % pad polar regions
    P{ii}=pad_sphgrid_poles(el,P{ii},configs.anal.maxPhi,NaN);
end


if dbg
    for ii=1:nparam
        this_parvec=parvec_tab{ii};   % get param-vec: [mf, traman]
        
        %%% plot
        figname=char(strjoin(string(this_parvec),'_'));
        figname=['P_',figname];
        
        figure('Name',figname);
        plotFlatMapWrappedRad(az,el,P{ii});
        cbar=colorbar('southoutside');
        cbar.Label.String='P';

        % save
        if exist('SAVEFIG','var') && SAVEFIG
            fname=fullfile(dir_base,'out','img',figname);
            print(fname,'-dpng');
        end
        
    end
end


%% 3. Rotation angle

% P_excited = sin(theta/2)^2
%   Theta   --> P
%   pi      --> 1
%   0       --> 0
%

th=cellfun(@(p) 2*asin(p),P,'UniformOutput',false);

if dbg
    for ii=1:nparam
        this_parvec=parvec_tab{ii};   % get param-vec: [mf, traman]
        
        %%% plot
        figname=char(strjoin(string(this_parvec),'_'));
        figname=['theta_',figname];
        
        figure('Name',figname);
        plotFlatMapWrappedRad(az,el,th{ii}/pi);
        cbar=colorbar('southoutside');
        cbar.Label.String='\theta/\pi';

        % save
        if exist('SAVEFIG','var') && SAVEFIG
            fname=fullfile(dir_base,'out','img',figname);
            print(fname,'-dpng');
        end
    end
end


%% 4. Clean data
% concatenate raman pulse duration
P_mf=cell(1,2);
for ii=1:2
    P_mf{ii}=cat(3,P{ii,:});
end

th_mf=cell(1,2);
for ii=1:2
    th_mf{ii}=cat(3,th{ii,:});
end


