%Analyse experiment on mF asymmetry
%

dbg=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BEGIN PREPROCESSING
addpath(genpath('C:\Users\HE BEC\exp-data\bell\loop_tscan'));


%% 1. Load raw data
% create txy and save mat file

% or LOAD IT!
load('txy_20180129.mat');
shot_id=fout.id_ok;


%% 1.1. Load parameter log
% see load_params.m
%   params
%   id_in_params

load('param_20180130.mat');

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
    this_par=params(ii,:);      % this param-vect
    this_idx=num2cell(par_tab(ii,:));     % location in table 
    %cell-formated to distribute as argument list
    
    k_par{this_idx{:}}=k_halo(b_paramset(:,ii),:);      % get all halo data
    %from this param-set and store
end

if dbg
    for ii=1:numel(k_par)
%         this_idx=cell(1,dim_param);
%         [this_idx{:}]=ind2sub(n_par_iter,ii);	% get array subscripts
        
        h=summary_disthalo_ndist(k_par{ii});
        
        
    end
end


%%% END OF PREPROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BEGIN ANALYSIS 

%% 1. 
