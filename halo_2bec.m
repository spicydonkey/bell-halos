function [halo_k,bec_cent]=halo_2bec(zxy,p_bec1,p_bec2,r_bec,r_th,dR_halo,zcap,verbose)
% captures halo based on 2 source BECs at poles
%
% DKS 31/10/2017
%
% TODO
% [] cull thermal fraction (r_th) around BEC center
% [] return more info
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
% capture BEC to mark halo poles
[bec_cent,bool_bec]=capture_bec(zxy,{p_bec1,p_bec2},{r_bec,r_bec},verbose);

% get atoms in BEC
% bool_bec_combined=cellfun(@(x,y)x|y,bool_bec(:,1),bool_bec(:,2),'UniformOutput',false);
bool_bec_combined=cell_horzcat(bool_bec);       % boolean array to indicate which BEC atoms belongs to
bool_bec_combined=cellfun(@(x)sum(x,2)>0,bool_bec_combined,'UniformOutput',false);    % row-wise OR
% boolean for: is this atom from a BEC?

%% 2. First stage halo capture
%%% 2.1. Locate halo center 
% get halo centre: midpoint of 'BEC poles'
halo_cent=cellfun(@(c1,c2)(c1+c2)/2,bec_cent(:,1),bec_cent(:,2),'UniformOutput',false);

% approximate halo radius: 0.5 * distance(BEC1,BEC2)
R_halo=cellfun(@(c1,c2) sqrt(sumsqr((c1-c2)/2)),bec_cent(:,1),bec_cent(:,2));
R_halo_mean=mean(R_halo);
R_halo_std=std(R_halo);

%%% 2.2. first capture of halo counts
% get counts in halo by popping BEC from raw
halo_zxy=cellfun(@(ZXY,BOOL)ZXY(~BOOL,:),zxy,bool_bec_combined,'UniformOutput',false);

% roughly center the halo to each shot-wise evaluated halo center
halo_zxy0=cellfun(@(zxy,c0)zxy-repmat(c0,[size(zxy,1),1]),halo_zxy,halo_cent,'UniformOutput',false);

%%% 2.3. filter halo
%%%% radial
% distance from approx halo centre
r0=cellfun(@(x)sqrt(sum(x.^2,2)),halo_zxy0,'UniformOutput',false);
rlim_hicap=R_halo_mean*(1+dR_halo*kdR_hicap*[-1,1]);       % radial limits for sph-shell capture

bool_halo_r_hicap=cellfun(@(r)(r<rlim_hicap(2))&(r>rlim_hicap(1)),r0,'UniformOutput',false);    % atoms in radial limits

% TODO: keep regional IN bool-array from each filtering stage, AND it when the
% filtered halo is needed
% % capture halo!
% halo_zxy0=cellfun(@(ZXY,BOOL)ZXY(BOOL,:),halo_zxy0,bool_halo_r_hicap,'UniformOutput',false);

%%%% elevation (TODO)
% removing caps using Z complicates angular analysis: use elevation angles
% instead
% % remove caps
% % TODO - try removing the caps after ellipsoid fit
% dz_poles=cellfun(@(c1,c2)abs(c1(1)-c2(1)),bec_cent(:,1),bec_cent(:,2),'UniformOutput',false);
% dz_poles=mean(vertcat(dz_poles{:}));    % pole-pole distance (BEC) in Z
% 
% bool_zcap=cellfun(@(zxy)abs(zxy(:,1))>(0.5*dz_poles*zcap),halo_zxy0,'UniformOutput',false);
% 
% halo_zxy0=cellfun(@(ZXY,BOOL)ZXY(~BOOL,:),halo_zxy0,bool_zcap,'UniformOutput',false);

%%%% thermal (TODO)
% centered around BEC locations but radius is ~2-3 times larger than used
% for BEC?

%%%% Filtered halo
% TODO - this is nasty joining cells like this. user must supply how many
% filters used: n_halo_filt1
% combine different filters
n_halo_filt1=1;
bool_halo_filt1_joined=cell(size(bool_halo_r_hicap,1),n_halo_filt1);
bool_halo_filt1_joined(:,1)=bool_halo_r_hicap;

% get first stage filter
bool_halo_filt1=cell_horzcat(bool_halo_filt1_joined);
bool_halo_filt1=cellfun(@(x) all(x,2),bool_halo_filt1,'UniformOutput',false);

% apply the filter
halo_zxy0=cellfun(@(ZXY,BOOL)ZXY(BOOL,:),halo_zxy0,bool_halo_filt1,'UniformOutput',false);

%% 3. Demanipulate halo
%%% 3.1. Ellipsoid fit
efit_flag='';
[ecent,erad,evecs,ev,echi2]=ellipsoid_fit(circshift(vertcat(halo_zxy0{:}),-1,2),efit_flag);

%%% 3.2. Unit sphere mapping 
%%%% Tranform to unit sphere (k-space)
% Initialise k-vector in XYZ coord
halo_k=cellfun(@(x) circshift(x,-1,2),halo_zxy0,'UniformOutput',false);

% Centre to fitted ellipsoid
halo_k(:)=boost_zxy(halo_k,-ecent');

% Transform to ellipsoid principal axis (Rotation)
M_rot=evecs;   % rotation matrix: X'Y'Z'(principal ellipsoid) --> XYZ
halo_k=cellfun(@(x) (M_rot\x')',halo_k,'UniformOutput',false);     % inverse transform

% Stretch to UNIT sphere: unit in collision wave-vector/momenta
halo_k=cellfun(@(x) x./repmat(erad',size(x,1),1),halo_k,'UniformOutput',false);

% Reverse transform to original/detector axis
halo_k=cellfun(@(x) (M_rot*x')',halo_k,'UniformOutput',false);

% transform to ZXY system
halo_k=cellfun(@(x) circshift(x,1,2),halo_k,'UniformOutput',false);


end