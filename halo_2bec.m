function [halo_k,bec_cent]=halo_2bec(zxy,p_bec1,p_bec2,r_bec,r_th,dR_halo,elev_max,verbose)
% captures halo based on 2 source BECs at poles
%
% DKS 31/10/2017
%
% TODO
% [x] cull thermal fraction (r_th) around BEC center
%   [x] r_th is an absolute magnitude now: need to be ~2*r_bec
%       [] update configs
%   [] test
% [] return more info
% [] comment code
%

% TODO
% [] improve algorithm
% 1. capture BEC: get BEC positions ==> get approx to halo center
% 2. 1st stage halo capture
%   2.1. liberal filtering of counts: radial, BEC, thermal
%   2.2. ellipsoid mapping --> unit sphere
% 3. clean halo
%   3.1. radial filter
%   3.2. pole filter
%

% - kdR_hicap: factor to make first stage radial culling more liberal


% parse input
if ~exist('verbose','var')
    warning('verbose is not provided - setting to quiet (0)');
    verbose=0;
end

%% Configs 
kdR_hicap=2;


%% 1. Capture marker BECs
%%% 1.1. capture BEC to mark halo poles
[bec_cent,bool_bec]=capture_bec(zxy,{p_bec1,p_bec2},{r_bec,r_bec},verbose);

%%% 1.2. get atoms in BEC
bool_bec_combined=cell_horzcat(bool_bec);       % boolean array to indicate which BEC atoms belongs to
bool_bec_combined=cellfun(@(b) any(b,2),bool_bec_combined,'UniformOutput',false);    % row-wise OR

%%% 1.3. thermal
% centered around BEC locations but radius can be ~2-3 times larger than used
% for BEC
bool_thermal=cell(size(zxy,1),2);
for ii=1:2
    % get BEC centered counts
    this_zxy0=cellfun(@(x,xc) x-repmat(xc,size(x,1),1),zxy,bec_cent(:,ii),'UniformOutput',false);
    this_r0=cellfun(@(x) zxy2rdist(x),this_zxy0,'UniformOutput',false);
    this_bool_thermal=cellfun(@(r) r<r_th,this_r0,'UniformOutput',false);
    this_bool_thermal=cellfun(@(b1,b2) b1&(~b2),this_bool_thermal,bool_bec_combined);      % select thermal but not BEC
    
    bool_thermal(:,ii)=this_bool_thermal;
end
clearvars this_zxy0 this_r0 this_bool_thermal;
bool_thermal_combined=cell_horzcat(bool_thermal);
bool_thermal_combined=cellfun(@(b) any(b,2),bool_thermal_combined,'UniformOutput',false);

%%% 1.4. collate BEC and thermal counts: a flare in image
% BEC and thermal counts are extremely large signals compared to scattered
% particle signal of interest
bool_flare=cellfun(@(b1,b2) b1|b2,bool_bec_combined,bool_thermal_combined,'UniformOutput',false);


%% 2. First stage halo capture
%%% 2.1. Locate halo center 
% get halo centre: midpoint of 'BEC poles'
halo_cent=cellfun(@(c1,c2)(c1+c2)/2,bec_cent(:,1),bec_cent(:,2),'UniformOutput',false);

% approximate halo radius: 0.5 * distance(BEC1,BEC2)
R_halo=cellfun(@(c1,c2) sqrt(sumsqr((c1-c2)/2)),bec_cent(:,1),bec_cent(:,2));
R_halo_mean=mean(R_halo);
R_halo_std=std(R_halo);

%%% 2.2. first capture of halo counts
% get counts in halo by popping BEC and thermal from raw
halo_zxy=cellfun(@(x,b) x(~b,:),zxy,bool_flare,'UniformOutput',false);

% roughly center the halo to each shot-wise evaluated halo center
halo_zxy0=cellfun(@(zxy,c0)zxy-repmat(c0,[size(zxy,1),1]),halo_zxy,halo_cent,'UniformOutput',false);

%%% 2.3. filter halo
%%%% radial
% distance from approx halo centre
r0=cellfun(@(x)sqrt(sum(x.^2,2)),halo_zxy0,'UniformOutput',false);
rlim_hicap=R_halo_mean*(1+dR_halo*kdR_hicap*[-1,1]);       % radial limits for sph-shell capture

bool_halo_r_hicap=cellfun(@(r)(r<rlim_hicap(2))&(r>rlim_hicap(1)),r0,'UniformOutput',false);    % atoms in radial limits

%%%% polar
% TODO - test does this work? cellfun multiple output?
[~,elev_halo_zxy0]=cellfun(@(zxy) zxy2sphpol(zxy),halo_zxy0,'UniformOutput',false);  % get elev angles [-pi/2,pi/2]
bool_halo_elev_hicap=cellfun(@(el) (abs(el)<elev_max),elev_halo_zxy0,'UniformOutput',false);    % atoms in elev limits

%%%% Filtered halo
% TODO - this is nasty joining cells like this. user must supply how many
% filters used: n_halo_filt1
% combine different filters
n_halo_filt1=2;
bool_halo_filt1_joined=cell(size(bool_halo_r_hicap,1),n_halo_filt1);
bool_halo_filt1_joined(:,1)=bool_halo_r_hicap;
bool_halo_filt1_joined(:,2)=bool_halo_elev_hicap;

% get first stage filter
bool_halo_filt1=cell_horzcat(bool_halo_filt1_joined);
bool_halo_filt1=cellfun(@(x) all(x,2),bool_halo_filt1,'UniformOutput',false);

% apply the filter
halo_zxy0=cellfun(@(z,b)z(b,:),halo_zxy0,bool_halo_filt1,'UniformOutput',false);

%% 3. Demanipulate halo
%%% 3.1. Ellipsoid fit
efit_flag='';
[ecent,erad,evec]=ellipsoid_fit(circshift(vertcat(halo_zxy0{:}),-1,2),efit_flag);   % collate all shots for fitting

%%% 3.2. Unit sphere mapping 
%%%% Tranform to unit sphere (k-space)
halo_k=cellfun(@(v_zxy) ellip2usph(v_zxy,ecent,erad,evec,verbose),...
    halo_zxy0,'UniformOutput',false);

%% 4. Clean the halo
%%% 4.1. Final filters
%%%% Radial
r_halo_k=cellfun(@(k_zxy) sqrt(sum(k_zxy.^2,2)),halo_k,'UniformOutput',false);
rlim_clean=1+dR_halo*[-1,1];        % a final hard crop
bool_halo_r_clean=cellfun(@(r) (r<rlim_clean(2))&(r>rlim_clean(1)),r_halo_k,'UniformOutput',false);    % atoms in radial limits

%%%% Polar
% TODO - check as above usage
[~,elev_halo_k]=cellfun(@(zxy) zxy2sphpol(zxy),halo_k,'UniformOutput',false);  % get elev angles [-pi/2,pi/2]
bool_halo_elev_clean=cellfun(@(el) (abs(el)<elev_max),elev_halo_k,'UniformOutput',false);    % atoms in elev limits

%%% 4.2. Filtered halo
% TODO - this is nasty joining cells like this. user must supply how many
% filters used: n_halo_clean
% combine different filters
n_halo_clean=2;
bool_halo_clean_joined=cell(size(bool_halo_r_clean,1),n_halo_clean);
bool_halo_clean_joined(:,1)=bool_halo_r_clean;
bool_halo_clean_joined(:,2)=bool_halo_elev_clean;

% get first stage filter
bool_halo_clean=cell_horzcat(bool_halo_clean_joined);
bool_halo_clean=cellfun(@(x) all(x,2),bool_halo_clean,'UniformOutput',false);

% apply the filter
halo_k=cellfun(@(z,b)z(b,:),halo_k,bool_halo_clean,'UniformOutput',false);


end