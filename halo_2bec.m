function [halo_k,bec_cent,hfig]=halo_2bec(zxy,p_bec1,p_bec2,r_bec,r_th,dR_halo,elev_max,verbose)
% captures halo based on 2 source BECs at poles
%
% [halo_k,bec_cent,hfig]=halo_2bec(zxy,p_bec1,p_bec2,r_bec,r_th,dR_halo,elev_max,verbose)
%
% DKS 31/10/2017
%
% 

% TODO
% [x] cull thermal fraction (r_th) around BEC center
%   [x] r_th is an absolute magnitude now: need to be ~2*r_bec
%       [x] update configs
%   [] test
% [] return more info
%   [] verbose output
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

% initialise figure outputs
hfig=[];    % figure array

%% 1. Capture marker BECs
%%% 1.1. capture BEC to mark halo poles
[bec_cent,bool_bec]=capture_bec(zxy,{p_bec1,p_bec2},{r_bec,r_bec},verbose-1);

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
    this_bool_thermal=cellfun(@(b1,b2) b1&(~b2),this_bool_thermal,bool_bec_combined,'UniformOutput',false);      % select thermal but not BEC
    
    bool_thermal(:,ii)=this_bool_thermal;
end
clearvars this_zxy0 this_r0 this_bool_thermal;
bool_thermal_combined=cell_horzcat(bool_thermal);
bool_thermal_combined=cellfun(@(b) any(b,2),bool_thermal_combined,'UniformOutput',false);

%%% 1.4. collate BEC and thermal counts: a flare in image
% BEC and thermal counts are extremely large signals compared to scattered
% particle signal of interest
bool_flare=cellfun(@(b1,b2) b1|b2,bool_bec_combined,bool_thermal_combined,'UniformOutput',false);

%%% Summary
% get all evaluated BEC cents
bec_cent_all=cell(1,2);
for ii=1:2
    bec_cent_all{ii}=vertcat(bec_cent{:,ii});
end
bec_cent_mean=cellfun(@(x) mean(x,1),bec_cent_all,'UniformOutput',false);
bec_cent_std=cellfun(@(x) std(x,1),bec_cent_all,'UniformOutput',false);

% output to screen
if verbose>0
    % raw zxy
    h_zxy_raw=figure();
    plot_zxy(zxy,1e4,1,'k');
    xlabel('X [m]');
    ylabel('Y [m]');
    zlabel('Z [m]');
    title('Raw ZXY');
    axis equal;
    box on;
    view(3);
end

if verbose>0
    % BEC cents
    for ii=1:2
        fprintf('[%s]: %s: bec_cent_mean{%d}=(%0.4g, %0.4g, %0.4g)\n',mfilename,getdatetimestr,ii,bec_cent_mean{ii});
        fprintf('[%s]: %s: bec_cent_std{%d}=(%0.4g, %0.4g, %0.4g)\n',mfilename,getdatetimestr,ii,bec_cent_std{ii});
%         fprintf('==========================================================================================\n');    %90 chars
    end
end
    
% oscillation compensated BEC/thermal categorized counts 
zxy0_flare=cellfun(@(x,b) x(b,:),zxy,bool_flare,'UniformOutput',false);
zxy0_remnant=cellfun(@(x,b) x(~b,:),zxy,bool_flare,'UniformOutput',false);
    
if verbose>0
    h_zxy0=figure();
    hold on;
    plot_zxy(zxy0_flare,1e4,1,'r');
    plot_zxy(zxy0_remnant,1e4,10,'k');
    xlabel('X [m]');
    ylabel('Y [m]');
    zlabel('Z [m]');
    title('Captured BEC/thermal');
    axis equal;
    box on;
    view(3);
    
    drawnow
end

% clean workspace
clearvars zxy0_flare zxy_remnant;


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
halo_zxy0_sphpol=cellfun(@(zxy) zxy2sphpol(zxy),halo_zxy0,'UniformOutput',false);  
% get elev angles [-pi/2,pi/2]
halo_zxy0_el=cellfun(@(pol) pol(:,2),halo_zxy0_sphpol,'UniformOutput',false);
bool_halo_elev_hicap=cellfun(@(el) (abs(el)<elev_max),halo_zxy0_el,'UniformOutput',false);    % atoms in elev limits
clearvars halo_zxy0_sphpol halo_zxy0_el;

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

%%% Summary
if verbose>0
    % Halo size
    fprintf('[%s]: %s: R_halo_mean=%0.3g\n',mfilename,getdatetimestr,R_halo_mean);
    fprintf('[%s]: %s: R_halo_std=%0.3g\n',mfilename,getdatetimestr,R_halo_std);
    %         fprintf('==========================================================================================\n');    %90 chars
end

if verbose>0
    % plot first stage captured halos
    h_halo_zxy0=figure();
    plot_zxy(halo_zxy0,1e4,10,'b');
    xlabel('X [m]');
    ylabel('Y [m]');
    zlabel('Z [m]');
    title('Halo (1st stage filtering)');
    axis equal;
    box on;
    view(3);
    drawnow
end

%% 3. Demanipulate halo
%%% 3.1. Ellipsoid fit
efit_flag='';
[ecent,erad,evec]=ellipsoid_fit(circshift(vertcat(halo_zxy0{:}),-1,2),efit_flag);   % collate all shots for fitting

%%% 3.2. Unit sphere mapping 
%%%% Tranform to unit sphere (k-space)
halo_k=cellfun(@(v_zxy) ellip2usph(v_zxy,ecent,erad,evec,verbose-1),...
    halo_zxy0,'UniformOutput',false);

%%% Summary
if verbose>0
    % unit sphere mapped halo
    h_halo_k=figure();
    plot_zxy(halo_k,1e4,10,'k');
    xlabel('$K_x$');
    ylabel('$K_y$');
    zlabel('$K_z$');
    title('Halo - unit-sphere mapped');
    axis equal;
    box on;
    view(3);
    drawnow;
end
    
%% 4. Clean the halo
%%% 4.1. Final filters
%%%% Radial
r_halo_k=cellfun(@(k_zxy) sqrt(sum(k_zxy.^2,2)),halo_k,'UniformOutput',false);
rlim_clean=1+dR_halo*[-1,1];        % a final hard crop
bool_halo_r_clean=cellfun(@(r) (r<rlim_clean(2))&(r>rlim_clean(1)),r_halo_k,'UniformOutput',false);    % atoms in radial limits

%%%% Polar
halo_k_sphpol=cellfun(@(zxy) zxy2sphpol(zxy),halo_k,'UniformOutput',false);
halo_k_el=cellfun(@(pol) pol(:,2),halo_k_sphpol,'UniformOutput',false);
bool_halo_elev_clean=cellfun(@(el) (abs(el)<elev_max),halo_k_el,'UniformOutput',false);    % atoms in elev limits
clearvars halo_k_sphpol halo_k_el;

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

%%% Summary
if verbose>0
    % draw over unit-sph mapped halo
    figure(h_halo_k); 
    hold on;
    %h_halo_k_filt=figure();
    plot_zxy(halo_k,1e4,20,'g');
    hold off;
    drawnow
end

%% 5. Radial distribution
% TODO
% - [x] radial distribution + Gaussian fit
%   - [] TEST!
%

% get radial dist to origin
halo_R=cellfun(@(c) zxy2rdist(c),halo_k,'UniformOutput',false);
halo_R=vertcat(halo_R{:});

%%% radial histogram
rhist_nbin=100;
%rhist_nbin=ceil(length(halo_R)/100);    % number of bins
[n_R,R_edge]=histcounts(halo_R,rhist_nbin);     % count in radial bins
R_cent=R_edge(1:end-1)+0.5*diff(R_edge);        % bin centers
n_R=n_R./(sum(n_R)*diff(R_edge));      % number to probability density function
% TODO - need to normalise 3D --> 1D radial (use R_cent)

if verbose>0
    h_rdist=figure();
    plot(R_cent,n_R,'o');   % radial histogram of data fitted to unit sphere

    box on;
    title('Halo radial distribution');
    xlabel('K');
    ylabel('PDF');

    drawnow
end

%%% fit distribution
% Model - 1D Gaussian with 0 offset
% see fit_gauss_1d

parameq={[],[],[],0};   % offset set to 0
param0=[10,1,0.033];   % [amplitude,mean,rmswidth]
[paramfit,~,rdist_gfit]=fit_gauss_1d(R_cent,n_R,param0,parameq);

% get fitted results
dk_rms=paramfit(3,:);   % rms halo thickness, fit SE
R_fit=linspace(min(R_cent),max(R_cent));
nR_fit=feval(rdist_gfit,R_fit);

if verbose>0
    fprintf('[%s]: %s: dk_rms=%0.3g (%0.1g)\n',mfilename,getdatetimestr,dk_rms);
end

if verbose>0
    figure(h_rdist);
    hold on;
    plot(R_fit,nR_fit,'k--','LineWidth',2);

    drawnow
end


%% 6. Number in halo
% TODO
% - [x] mean and sdev 
% - [] account for filters
%   - [] elev filter
%   - [x] QE?
% - [] TEST

ncounts_halo_k_shots=cellfun(@(k) size(k,1),halo_k);  % counts in halo_k per shot
n_scat=[mean(ncounts_halo_k_shots),std(ncounts_halo_k_shots)];  % number in filtered halo

% estimate total number in scattered halo (unfiltered)
vol_factor_elev=1;      % TODO - formula for this? f(elev_max)
det_qe=0.1;     % TODO may be closer to 0.08?
n_scat=n_scat/(det_qe*vol_factor_elev);

%%% Summary
if verbose>0
    fprintf('[%s]: %s: n_scat=%0.3g (%0.1g)\n',mfilename,getdatetimestr,n_scat);
end


%% Collate figures
hfig=[h_zxy_raw,h_zxy0,h_halo_zxy0,h_halo_k,h_rdist];

end
