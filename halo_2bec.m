function [halo_k,bec_cent]=halo_2bec(zxy,p_bec1,p_bec2,r_bec,r_th,dR_halo,zcap,verbose)
% captures halo based on 2 source BECs at poles
%
% DKS 31/10/2017
%
% TODO
% [] cull thermal fraction (r_th) around BEC center
% [] return more info
%

% parse input
if ~exist('verbose','var')
    warning('verbose is not provided - setting to quiet (0)');
    verbose=0;
end


%% Capture BECs and locate halo center
% capture BEC to mark halo poles
[bec_cent,bool_bec]=capture_bec(zxy,{p_bec1,p_bec2},{r_bec,r_bec},verbose);

% pop bec + thermal counts
% bool_bec_combined=cellfun(@(x,y)x|y,bool_bec(:,1),bool_bec(:,2),'UniformOutput',false);
bool_bec_combined=cell_horzcat(bool_bec);       % horizontal concat the boolean list
bool_bec_combined=cellfun(@(x)sum(x,2)>0,bool_bec_combined,'UniformOutput',false);    % row-wise OR

halo_zxy=cellfun(@(ZXY,BOOL)ZXY(~BOOL,:),zxy,bool_bec_combined,'UniformOutput',false);

% halo centre: midpoint of 'BEC poles'
halo_cent=cellfun(@(c1,c2)(c1+c2)/2,bec_cent(:,1),bec_cent(:,2),'UniformOutput',false);
% halo radius
R_halo=cellfun(@(c1,c2) sqrt(sumsqr((c1-c2)/2)),bec_cent(:,1),bec_cent(:,2));
R_halo_mean=mean(R_halo);
R_halo_std=std(R_halo);

% centre to halo
halo_zxy0=cellfun(@(zxy,c0)zxy-repmat(c0,[size(zxy,1),1]),halo_zxy,halo_cent,'UniformOutput',false);


%% radial capture
% distance from halo centre
r0=cellfun(@(x)sqrt(sum(x.^2,2)),halo_zxy0,'UniformOutput',false);
rlim=R_halo_mean*(1+dR_halo*[-1,1]);       % radial limits for sph-shell capture

bool_halo=cellfun(@(r)(r<rlim(2))&(r>rlim(1)),r0,'UniformOutput',false);    % determine halo points

% capture halo!
halo_zxy0=cellfun(@(ZXY,BOOL)ZXY(BOOL,:),halo_zxy0,bool_halo,'UniformOutput',false);

% remove caps
% TODO - try removing the caps after ellipsoid fit
dz_poles=cellfun(@(c1,c2)abs(c1(1)-c2(1)),bec_cent(:,1),bec_cent(:,2),'UniformOutput',false);
dz_poles=mean(vertcat(dz_poles{:}));    % pole-pole distance (BEC) in Z

bool_zcap=cellfun(@(zxy)abs(zxy(:,1))>(0.5*dz_poles*zcap),halo_zxy0,'UniformOutput',false);

halo_zxy0=cellfun(@(ZXY,BOOL)ZXY(~BOOL,:),halo_zxy0,bool_zcap,'UniformOutput',false);


%% Elipsoid fit
efit_flag='';
[ecent,erad,evecs,ev,echi2]=ellipsoid_fit(circshift(vertcat(halo_zxy0{:}),-1,2),efit_flag);


%% Tranform to unit sphere (k-space)
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